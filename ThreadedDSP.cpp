#include "ThreadedDSP.hpp"
#include <stdlib.h>
#include <assert.h>
#include <cstring>
#include <math.h>

#ifdef __linux__
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/fcntl.h>  
#include <sys/shm.h> 
#define mbstowcs_s(a,b,c,d,e) mbstowcs(b,d,e);
#else
#ifndef UNICODE
#define UNICODE
#endif
#include <Windows.h>
#include <AVRt.h>
#undef min
#undef max
#endif

ThreadedDSP::ThreadedDSP()
{
    hMem = 0;
    p    = 0;
    nSize = 0;
    bOwner   = false;
	pFFTSpec = pFFTInit = pFFTBuf  = 0;
	pFFT = 0;
    LastT = 0;
    done = 0;
    processing = false;
    pfX = 0;
    pfY = 0;
}

bool ThreadedDSP::Create(int block, int blocks, int in, int out, int filters, int xfade, const char* sShared, int threads, int latency)
{
    if (hMem || p)  return false;           // We are already created or attached

    if (block<32  || block>MaxBlock)           return false;
    if (blocks<4  || blocks>(MaxLength/block)) return false;
    if (in<0      || in>MaxChannels)           return false;
    if (out<0     || out>MaxChannels)          return false;
    if (filters<0 || filters>MaxFilters)       return false;
    if (sShared == nullptr)                    return false;
    if (threads<1 || threads>32)               return false;

    int N = block + xfade / 2;
    N--; N |= N >> 1; N |= N >> 2; N |= N >> 4; N |= N >> 8; N++;   // Make it the next power of 2
    xfade = std::min(2*(N - block),block);

    printf("CREATING CONVOLVER N=%d  BLOCK=%d  XFADE=%d    BURDEN=%.0f%%   UTILIZATION=%.0f%%\n",N,block,xfade,100.0F-100.0F*block*block/N/N,100.0F*(block+block+xfade)/(2*N));

    int FFTorder = (int)(logf((float)N) / logf(2.0F) + 1.5);
  
#ifdef __linux__
    if (sShared[0]!='/') { sMem[0]='/'; strncpy(sMem+1,sShared,254); } else strncpy(sMem, sShared, 255);
    mode_t old_mask = umask(0);
    hMem = (void *)(long)shm_open(sMem, O_RDWR | O_CREAT | O_EXCL,  S_IRWXU | S_IRWXG | S_IRWXO);
    umask(old_mask);
    if (hMem==(void *)(long)(-1)) { hMem = 0; sMem[0]=0; return false; };
    nSize = Size(in, out, blocks, block+xfade/2, filters, block);
    if (ftruncate((int)(long)hMem,nSize)==-1) { shm_unlink(sMem); hMem=0; sMem[0]=0; nSize=0; return false; };
    p = (ThreadedDSP::Map *)mmap(NULL,nSize,PROT_READ|PROT_WRITE,MAP_SHARED, (int)(long)hMem, 0);
    if (p==NULL) { shm_unlink(sMem); hMem=0; sMem[0]=0; nSize=0; return false; };
#else
    wchar_t wc[64];
    mbstowcs_s(0, wc, 63, sShared, 63);
    hMem = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, 0, wc);
    if (hMem)                                         {  CloseHandle(hMem);   assert(0); };
    nSize = Size(in, out, blocks, block+xfade/2, filters);
    hMem = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, (DWORD)nSize, wc);
    if (hMem && GetLastError()==ERROR_ALREADY_EXISTS) {  CloseHandle(hMem);   hMem=0; nSize=0; return false; };
    if (hMem==NULL) { nSize=0; return false; };
    p = (Map *)MapViewOfFile(hMem, FILE_MAP_ALL_ACCESS, 0, 0, nSize);
    if (p==NULL) { CloseHandle(hMem); hMem=0; nSize=0; return false; };
#endif

    memset(p,0,nSize);
    bOwner = true;
    p->Block = block; 
    p->N     = N; 
    p->M     = blocks; 
    p->I     = in; 
    p->O     = out;
    p->F     = filters;
    p->MBlock= p->M*p->Block;
    p->MN    = p->M*p->N;
    p->Latency = latency;
    p->Threads = threads;
    nThreads = threads;     // Keep a local copy for the owner so it cannot be changed

    pfX = (float *)malloc(2*p->I*p->M*p->N*sizeof(float));
    pfY = (float *)malloc(2*p->O*p->N*sizeof(float));

    pW  = (float *)malloc((2*p->N - p->Block)*sizeof(float));

    for (int n=0;        n<xfade;    n++) pW[n]          = 0.5F - 0.5F*cosf(3.14159265358979323846F*((float)n+0.5)/xfade);
    for (int n=xfade;    n<p->Block; n++) pW[n]          = 1.0F;
    for (int n=0;        n<xfade;    n++) pW[n+p->Block] = 0.5F + 0.5F*cosf(3.14159265358979323846F*((float)n+0.5)/xfade);

    for (int n=0; n<p->I; n++) p->GainIn[n] = 1.0F;
    for (int n=0; n<p->O; n++) p->GainOut[n] = 1.0F;

    //for (int n=0; n<2*N-p->Block; n++) printf("%2d %5.3f\n",n,pW[n]);
    
    anTemp = (int32_t *)malloc(p->N*sizeof(int32_t));

	int sizeSpec, sizeInit, sizeBuf;
    ippsFFTGetSize_R_32f(FFTorder, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, &sizeSpec, &sizeInit, &sizeBuf);
	pFFTSpec = (Ipp8u *)ippsMalloc_8u(sizeSpec);
	pFFTInit = (Ipp8u *)ippsMalloc_8u(sizeInit);
	pFFTBuf  = (Ipp8u *)ippsMalloc_8u(sizeBuf);
	ippsFFTInit_R_32f(&pFFT, FFTorder, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, pFFTSpec, pFFTInit);
    p->bRunning = true;

    float range = 3*block/48000.0F;

    Chrono_CallTime().configure("DSP CALL TIME",0.000F, range);
    Chrono_Load().configure("DSP EXECUTION LOAD %",0, 200);

    done = 0;
    for (int t=0; t<nThreads; t++) 
    { 
        char name[64];
        snprintf(name, 63, "DSP EXECUTION TIME  THREAD %2d", t);
        Chrono_ExecTime(t).configure(name, 0.000F, range);
        pThreads[t] = new std::thread(&ThreadedDSP::ThreadProc,this,t); 
    };
    std::this_thread::sleep_for(std::chrono::milliseconds(10)); // Make sure all threads start and get to first lock
    return p->bRunning;
}

bool ThreadedDSP::Attach(const char* sShared)
{
    if (hMem || p)  return false;           // We are already created or attached

#ifdef __linux__
    if (sShared[0]!='/') { sMem[0]='/'; strncpy(sMem+1,sShared,254); } else strncpy(sMem, sShared, 255);
    mode_t old_mask = umask(0);
    hMem = (void *)(long)shm_open(sMem, O_RDWR,            S_IRUSR | S_IWUSR);
    umask(old_mask);
    if (hMem==(void *)(long)(-1)) {  sMem[0]=0; return false; };      // Does not exist
    nSize = lseek((int)(long)hMem, 0, SEEK_END);
    p = (ThreadedDSP::Map *)mmap(NULL,nSize,PROT_READ|PROT_WRITE,MAP_SHARED, (int)(long)hMem, 0);
    if (p==NULL)                                               { shm_unlink(sMem); hMem=0; sMem[0]=0; nSize=0; return false; };
    if (nSize != Size(p->I, p->O, p->M, p->N, p->F, p->Block)) { shm_unlink(sMem); hMem=0; sMem[0]=0; nSize=0; return false; };
#else
    wchar_t wc[64];
    mbstowcs_s(0, wc, 63, sShared, 63);
    hMem = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, sizeof(Map), wc);       // Just open the length of the Map structure
    if (hMem==0) return false;                                                                      // to see if it exists
    p = (Map *)MapViewOfFile(hMem, FILE_MAP_ALL_ACCESS, 0, 0, 0);                                   //
    if (p==NULL || p->bRunning==false) { CloseHandle(hMem); return false; };                        //
    nSize = Size(p->I, p->O, p->M, p->N, p->F, p->Block);                                           // Get the full size
    UnmapViewOfFile(p); CloseHandle(hMem);                                                          // Close this
    hMem = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, (DWORD)nSize, wc);      // Then open it at the full size
    if (hMem==NULL) { nSize=0; return false; };
    p = (Map *)MapViewOfFile(hMem, FILE_MAP_ALL_ACCESS, 0, 0, nSize);
    if (p==NULL) { CloseHandle(hMem); hMem=0; nSize=0; return false; };
#endif
    nThreads = p->Threads;
    anTemp = (int32_t *)malloc(  p->N*sizeof(int32_t));
    return true;
}


ThreadedDSP::~ThreadedDSP(void)
{
    if (p) Destroy();
}

void ThreadedDSP::Destroy(void)
{
    if (hMem==nullptr) return;
    if (bOwner)
    {
        if (p) p->bRunning = false;
        start.notify_all(); std::this_thread::sleep_for(std::chrono::milliseconds(10));
        start.notify_all(); std::this_thread::sleep_for(std::chrono::milliseconds(10));

        for (int g=0; g<MaxGroups; g++) 
        {
            for (int s=0; s<MaxSets; s++) 
            {
                if (Filt_Set[g][s]) 
                {
                    if (s!=p->Filt_Active[g])
                    {
                        for (int f=0; f<p->F; f++)
                        {
                            if (Filt_Set[g][s][f].H) { free(Filt_Set[g][s][f].H);  Filt_Set[g][s][f].H = 0; };
                        }
                    }
                    delete [] Filt_Set[g][s];
                    Filt_Set[g][s] = 0;
                }
            }
            for (int f=0; f<p->F; f++)  if (p->Filt[g][f].H) { free(p->Filt[g][f].H);  p->Filt[g][f].H = 0; };
        }
    }
#ifdef __linux__
    if (p && nSize>0)                  munmap(p, nSize);
    if (bOwner && hMem && sMem[0] !=0) shm_unlink(sMem);
    p       = nullptr;
    hMem    = nullptr;
    sMem[0] = 0;
    nSize   = 0;
#else
    if (p)  	UnmapViewOfFile(p);     p        = nullptr;
	if (hMem)	CloseHandle(hMem);      hMem     = nullptr;
#endif
	if (pFFTBuf)  ippsFree(pFFTBuf);    pFFTBuf  = nullptr;
	if (pFFTInit) ippsFree(pFFTInit);   pFFTInit = nullptr;
    if (pFFTSpec) ippsFree(pFFTSpec);   pFFTSpec = nullptr;
    if (anTemp)   free(anTemp);         anTemp   = nullptr;

    if (pfX)      free(pfX);            pfX      = nullptr;
    if (pfY)      free(pfY);            pfY      = nullptr;
    if (pW)       free(pW);             pW       = nullptr;

    bOwner = false;

}


void ThreadedDSP::Input(int i, int32_t* data, int stride)
{
    if (!bOwner || i<0 || i>=p->I || !p || !p->bRunning) return;
    std::unique_lock<std::mutex> lk(mtx);
    if (!loading) { loading=true; p->T++; p->t = (p->t+1)%p->M; assert(p->T % p->M == p->t); }; // Advance
    lk.unlock();
    float   *pfIn   = afIn(i, 2*p->N - p->Block);
    
    if (p->GainIn[i] == 1.0F) for (int n=0; n<p->Block; n++)  {  *pfIn  = (float)*data / 2147483648.0F; pfIn++; data+=stride;  }
    else                      for (int n=0; n<p->Block; n++)  {  *pfIn  = (float)*data / 2147483648.0F * p->GainIn[i];  pfIn++; data+=stride; }
    // TODO The gain change here on input may cause zippering )
}

void ThreadedDSP::Output(int o, int32_t *data, int stride) 
{ 
    if (!bOwner || o<0 || o>=p->O || !p || !p->bRunning) return;
    if (data == 0)
    {
        return;
    }

    int32_t *pN  = anTemp;
    ippsConvert_32f32s_Sfs(afOut(o,p->Block), pN, p->Block, ippRndZero, -31);
    for (int n=0; n<p->Block; n++) 
    {
        *data = *pN++;
        data += stride;
    }
}

void ThreadedDSP::Input(int n, float *data, int stride)  { assert(0); }
void ThreadedDSP::Output(int n, float *data, int stride) { assert(0); }


void ThreadedDSP::Finish(void)
{
    if (!p || !bOwner) return;
    std::unique_lock<std::mutex> lk(mtx);           // Wait for all DSP threads to finish
    if (processing && p->bRunning) finished.wait(lk); lk.unlock(); // Only wait if not already finished
    float tmax = 0;
    for (int n=0; n<nThreads; n++) tmax = std::max(tmax, Chrono_ExecTime(n).latest());
    p->Load = (float)tmax / p->Block * nRate;
    Chrono_Load().add((int)(p->Load*100+0.5));
}

void ThreadedDSP::Process(void)
{
    if (!p || !bOwner) return;
    Chrono_CallTime().time();

    if (p->ClearAll) 
    {
        for (int g=0; g<MaxGroups; g++) 
        {
            for (int s=0; s<MaxSets; s++) 
            {
                if (Filt_Set[g][s]) 
                {
                    if (s!=p->Filt_Active[g])
                    {
                        for (int f=0; f<p->F; f++)
                        {
                            if (Filt_Set[g][s][f].H) { free(Filt_Set[g][s][f].H);  Filt_Set[g][s][f].H = 0; };
                        }
                    }
                    delete [] Filt_Set[g][s];
                    Filt_Set[g][s] = 0;
                }
            }
            for (int f=0; f<p->F; f++)  if (p->Filt[g][f].H) { free(p->Filt[g][f].H);  p->Filt[g][f].H = 0; p->Filt[g][f].BLen = 0; };
        }
        p->ClearAll = false;
    }
    
    for (int g=0; g<MaxGroups; g++)
    {
        if (p->Filt_Next[g] != p->Filt_Active[g]) 
        {
            if (Filt_Set[g][p->Filt_Active[g]]==0) Filt_Set[g][p->Filt_Active[g]] = new Filter[p->F]();
            if (Filt_Set[g][p->Filt_Next[g]]==0)   Filt_Set[g][p->Filt_Next[g]]   = new Filter[p->F]();
            for (int f=0; f<p->F; f++) { p->Filt[g][f].Update = false; Filt_Set[g][p->Filt_Next[g]][f].Update = false; };    // Avoid stale or lagging updates
            memcpy(Filt_Set[g][p->Filt_Active[g]], p->Filt, p->F*sizeof(Filter));
            memcpy(p->Filt, Filt_Set[g][p->Filt_Next[g]], p->F*sizeof(Filter));
            p->Filt_Active[g] = p->Filt_Next[g];
        }
    }
    for (int n=0; n<nThreads; n++) Chrono_ExecTime(n).start();  
    std::unique_lock<std::mutex> lk(mtx); loading=false; processing=true; start.notify_all(); lk.unlock();  // Start the DSP threads
}

void ThreadedDSP::ThreadProc(int ID)
{
#ifdef __linux__
	struct sched_param param{};
	param.sched_priority = 60;
    int ret = pthread_setschedparam(pthread_self(), SCHED_FIFO, &param);
#else
    DWORD taskIndex = 0;
	HANDLE h = AvSetMmThreadCharacteristics(TEXT("Pro Audio"), &taskIndex);
    SetPriorityClass(GetCurrentProcess(),HIGH_PRIORITY_CLASS);
    SetThreadPriority(GetCurrentThread(),THREAD_PRIORITY_HIGHEST);
    timeBeginPeriod(1);
#endif
    int sizeSpec, sizeInit, sizeBuf;
	float *afTemp = (float *)malloc(2*p->N*sizeof(float));
    float fMax;
	ippSetDenormAreZeros(true);
    int FFTorder = (int)(logf((float)p->N) / logf(2.0F) + 1.5);
	ippsFFTGetSize_R_32f(FFTorder, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, &sizeSpec, &sizeInit, &sizeBuf);
	Ipp8u *pFFTSpec = (Ipp8u *)ippsMalloc_8u(sizeSpec);
	Ipp8u *pFFTInit = (Ipp8u *)ippsMalloc_8u(sizeInit);
	Ipp8u *pFFTBuf  = (Ipp8u *)ippsMalloc_8u(sizeBuf);
	IppsFFTSpec_R_32f *pFFT;
	ippsFFTInit_R_32f(&pFFT, FFTorder, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, pFFTSpec, pFFTInit);


    std::unique_lock<std::mutex> lk (mtx);                              // Loop always starts locked
    while(p->bRunning)
    {
        start.wait(lk); lk.unlock();                                    // Wait on CV signal to start processing
        for (int n=ID; n<p->I; n+=nThreads) 			                // Do this over all of the inputs
		{				
            ippsMaxAbs_32f(afIn(n,2*p->N-p->Block),p->Block, &fMax);               // Get the peak
            ippsFFTFwd_RToPerm_32f(afIn(n,0), afX(n,p->t,0), pFFT, pFFTBuf);        // Do the FFT
            ippsMove_32f(afIn(n,p->Block), afIn(n,0), 2*p->N-p->Block);             // Slide the input data buffer back
            if (fMax>p->PeakIn[n]) p->PeakIn[n]=fMax; else p->PeakIn[n] = 0.99F*p->PeakIn[n];  // Bit of smoothing
        }
        lk.lock(); if (++done==nThreads) middle.notify_all(); else middle.wait(lk); lk.unlock();    // Resync

        if (!p->bRunning) break;
        for (int n=ID; n<p->O; n+=nThreads) 			                // Do this over all of the outputs
		{				
			ippsSet_32f (0.0F, afY(n,0), 2*p->N);  				        // Zero MAC buffer
			for (int g=0; g<MaxGroups; g++)
            {
                for (int f=0; f<p->F; f++) 
                {														//
                    if (p->Filt[g][f].Out != n) continue;		   		// Only filters with output channel n
                    if (p->Filt[g][f].Update) UpdateFilter(f, afTemp, pFFT, pFFTBuf, g);             // Update the filter - out wont change, but m may
                    for (int m=0; m<p->Filt[g][f].BLen; m++) 	        // Loop over the filter M blocks
                    {
                        int nB = (p->t - m + p->M)%p->M;       		    // and do the MAC
                        float r0=*afY(n,0), r1=*afY(n,1);               // Keep this as a multiple of 4 for SIMD optimization
                        ippsAddProduct_32fc((const Ipp32fc *)afX(p->Filt[g][f].In,nB,0),(const Ipp32fc *)afH(f,m,0), (Ipp32fc *)afY(n,0), p->N);
                        *afY(n,0) = r0 + *afX(p->Filt[g][f].In,nB,0) * *afH(f,m,0);    // Replace the first bin (DC + Nyquist)
                        *afY(n,1) = r1 + *afX(p->Filt[g][f].In,nB,1) * *afH(f,m,1);    // with the normal multiply
                    }
                }
            }
 //         ippsMove_32f(afX(n,p->t,0), afY(n,0), 2*p->N);                  // FFT LOOPBACK
            ippsMove_32f(afOut(n,p->Block), afOut(n,0), 2*p->N-p->Block);	// Slide the output data buffer back
            ippsFFTInv_PermToR_32f(afY(n,0), afTemp, pFFT,pFFTBuf);		    // Do the inverse FFT and stick it in the
            if (p->GainOut[n] != 1.0F) ippsMulC_32f_I(p->GainOut[n], afTemp+p->Block, 2*p->N-p->Block); // Apply the gain
            if (p->N != p->Block)
            {
             	ippsMul_32f_I(pW, afTemp+p->Block, 2*p->N-p->Block);		            // Apply the window
                ippsAdd_32f_I(afTemp+p->Block,afOut(n,p->Block),2*p->N-2*p->Block);     // Add the overlap
            }
            ippsMove_32f(afTemp+2*p->N-p->Block, afOut(n,2*p->N-p->Block),p->Block);    // Move the new data and fade out bit
//          ippsMove_32f(afIn(n,p->Block),afOut(n,p->Block),p->Block);      // INPUT AUDIO LOOPBACK
            ippsMaxAbs_32f(afOut(n,p->Block),p->Block,&fMax);               // Get the peak
            if (fMax>p->PeakOut[n]) p->PeakOut[n]=fMax; else p->PeakOut[n] = 0.99F*p->PeakOut[n];  // Bit of smoothing
		}
        Chrono_ExecTime(ID).time();
        lk.lock(); if (--done==0) { finished.notify_one(); processing=false; };// Check if this is the last, and return to start without unlocking
    }
    ippsFree(pFFTBuf); 
	ippsFree(pFFTInit);
    ippsFree(pFFTSpec);
    free(afTemp);
 }

void ThreadedDSP::UpdateFilter(int f,  float* afTemp, IppsFFTSpec_R_32f *pFFT, Ipp8u *pFFTBuf, int g)
{
    if (!bOwner || !p->Filt[g][f].Update) return;
    int length   = p->Filt[g][f].Length;

    int M = (length-1)/p->Block + 1;
    if (length == 0) M = 0;

    if (M != p->Filt[g][f].BLen) 
    {
        if (p->Filt[g][f].H) free(p->Filt[g][f].H);
        p->Filt[g][f].H = 0;
        if (M>0) p->Filt[g][f].H = (float *)calloc(p->N*M,2*sizeof(float));
    }

    float* pFilt = afT(f,0);
    int        m = 0;
	while (length > 0 && m < p->M)
	{
		float *pF = afTemp;
		int     n = p->Block;                   // Block coefficients in each filter block
		ippsSet_32f(0, afTemp, 2*p->N);         // But 2*N read samples for the FFT updload
		while (n > 0 && length > 0) { *pF++ = *pFilt++; length--; n--; };
 		ippsFFTFwd_RToPerm_32f(afTemp, afH(f,m,0), pFFT, pFFTBuf);
		m++;
    }
	p->Filt[g][f].BLen = m;
    p->Filt[g][f].Update = false;
}


bool ThreadedDSP::LoadFilter(int in, int out, int length, float *pFilt, int g)
{
    if (!p || !p->bRunning) return false;        
	if (in < 0 || in >= p->I || out < 0 || out >= p->O) return false;	
    
    int f, m=0;
    for (f = 0; f<p->F; f++) 
        if (p->Filt[g][f].In == in && p->Filt[g][f].Out == out) break;  // Check if we have this in->out already
    if (f == p->F)                                                      // 
    {                                                                   // If we don't have it
        if (length == 0) return true;                                   //   And if it is being cleared this is OK                                
    	for (f = 0; f<p->F; f++)                                        //
            if (p->Filt[g][f].Length==0 && p->Filt[g][f].BLen==0) break;//   Otherwise search for a new slot
        if (f == p->F) return false;                                    //   If none available then return false
    }                                                                   // At this point we have an existing or new slot
    if (length>p->MBlock) length = p->MBlock;
    if (length>0) memcpy(afT(f,0), pFilt, length*sizeof(float));
    assert(p->Filt[g][f].BLen == 0 || (p->Filt[g][f].In==in && p->Filt[g][f].Out==out));   // Just to be sure - not swapping live in/out
    p->Filt[g][f].In     = in;                 
    p->Filt[g][f].Out    = out;
    p->Filt[g][f].Length = length;
    p->Filt[g][f].Update = true;
    return true;
}

uint64_t ThreadedDSP::Spin(uint64_t T, int timeoutms, bool spin)
{
    if (!p || !p->bRunning || T<0 || p->T==0 || timeoutms<0 || timeoutms>10000) return 0;
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    if      (T==0)            T = p->T  + 1;
    if      (T<=p->M)         T = LastT + T;
    else if (T > p->T + p->M) T = p->T  + p->M;
    int t = *(volatile uint64_t *)&p->T;
    while ((T > t) && p->bRunning)
    {
        if (timeoutms>0 
            && std::chrono::duration_cast<std::chrono::milliseconds>
               (std::chrono::high_resolution_clock::now() - start).count() > timeoutms ) break;
        if (!spin && (T-*(volatile uint64_t *)&p->T)>1) std::this_thread::sleep_for(std::chrono::microseconds(100));
        // Note this will spin if there is less than one frame to wait - as a sleep is usually more than one frame
        t = *(volatile uint64_t *)&p->T;
    }
    LastT = T;
    return  LastT;
}
