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
}

bool ThreadedDSP::Create(int block, int blocks, int in, int out, int filters, int xfade, const char* sShared, int threads, int latency)
{
    assert(xfade==0);                       // Not yet supported
    if (hMem || p)  return false;           // We are already created or attached

    if (block<32  || block>MaxBlock)           return false;
    if (blocks<4  || blocks>(MaxLength/block)) return false;
    if (in<0      || in>MaxChannels)           return false;
    if (out<0     || out>MaxChannels)          return false;
    if (filters<0 || filters>MaxFilters)       return false;
    if (sShared == nullptr)                    return false;
    if (threads<1 || threads>32)               return false;

    int N = block+xfade/2;
   // assert((logf((float)N)/logf(2.0F))-(int)(logf((float)N)/logf(2.0F)+0.5));       // FFT block size must be power of 2 for IPP
    FFTorder = (int)(logf((float)N)/logf(2.0F)+1.5);
  
#ifdef __linux__
    if (sShared[0]!='/') { sMem[0]='/'; strncpy(sMem+1,sShared,254); } else strncpy(sMem, sShared, 255);
    mode_t old_mask = umask(0);
    hMem = (void *)(long)shm_open(sMem, O_RDWR | O_CREAT | O_EXCL,  S_IRWXU | S_IRWXG | S_IRWXO);
    umask(old_mask);
    if (hMem==(void *)(long)(-1)) { hMem = 0; sMem[0]=0; return false; };
    nSize = Size(in, out, blocks, block+xfade/2, filters);
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
    p->N     = block+xfade/2; 
    p->M     = blocks; 
    p->I     = in; 
    p->O     = out;
    p->F     = filters;
    p->MBlock= p->M*p->Block;
    p->MN    = p->M*p->N;
    p->Latency = latency;
    p->Threads = threads;
    nThreads = threads;     // Keep a local copy for the owner so it cannot be changed
    afTemp = (float *)malloc(2*p->N*sizeof(float));
    anTemp = (int32_t *)malloc(p->N*sizeof(int32_t));

    Filt_Set[0] = new Filter[p->F]();

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
    if (p==NULL)                                     { shm_unlink(sMem); hMem=0; sMem[0]=0; nSize=0; return false; };
    if (nSize != Size(p->I, p->O, p->M, p->N, p->F)) { shm_unlink(sMem); hMem=0; sMem[0]=0; nSize=0; return false; };
#else
    wchar_t wc[64];
    mbstowcs_s(0, wc, 63, sShared, 63);
    hMem = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, sizeof(Map), wc);       // Just open the length of the Map structure
    if (hMem==0) return false;                                                                      // to see if it exists
    p = (Map *)MapViewOfFile(hMem, FILE_MAP_ALL_ACCESS, 0, 0, 0);                                   //
    if (p==NULL || p->bRunning==false) { CloseHandle(hMem); return false; };                        //
    nSize = Size(p->I, p->O, p->M, p->N, p->F);
    UnmapViewOfFile(p); CloseHandle(hMem);                                                          // Close this
    hMem = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, (DWORD)nSize, wc);      // Then open it at the full size
    if (hMem==NULL) { nSize=0; return false; };
    p = (Map *)MapViewOfFile(hMem, FILE_MAP_ALL_ACCESS, 0, 0, nSize);
    if (p==NULL) { CloseHandle(hMem); hMem=0; nSize=0; return false; };
#endif
    nThreads = p->Threads;
    afTemp = (float   *)malloc(2*p->N*sizeof(float));
    anTemp = (int32_t *)malloc(  p->N*sizeof(int32_t));
	int sizeSpec, sizeInit, sizeBuf;
    FFTorder = (int)(logf((float)p->N)/logf(2.0F)+1.5);
    ippsFFTGetSize_R_32f(FFTorder, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, &sizeSpec, &sizeInit, &sizeBuf);
	pFFTSpec = (Ipp8u *)ippsMalloc_8u(sizeSpec);
	pFFTInit = (Ipp8u *)ippsMalloc_8u(sizeInit);
	pFFTBuf  = (Ipp8u *)ippsMalloc_8u(sizeBuf);
	ippsFFTInit_R_32f(&pFFT, FFTorder, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, pFFTSpec, pFFTInit);
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

        for (int s=0; s<MaxSets; s++) 
        {
            if (s==p->Filt_Active) continue;
            if (Filt_Set[s]) 
            {
                for (int f=0; f<p->F; f++)  if (Filt_Set[s][f].H) { free(Filt_Set[s][f].H);  Filt_Set[s][f].H = 0; };
                delete [] Filt_Set[s];
            }
        }
        for (int f=0; f<p->F; f++)  if (p->Filt[f].H) { free(p->Filt[f].H);  p->Filt[f].H = 0; };
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
    if (afTemp)   free(afTemp);         afTemp   = nullptr;
    if (anTemp)   free(anTemp);         anTemp   = nullptr;
    bOwner = false;

}


void ThreadedDSP::Input(int i, int32_t* data, int stride)
{
    if (!bOwner || i<0 || i>=p->I || !p || !p->bRunning) return;
    std::unique_lock<std::mutex> lk(mtx);
    if (!loading) { loading=true; p->T++; p->t = (p->t+1)%p->M; assert(p->T % p->M == p->t); }; // Advance
    lk.unlock();
    float   *pfIn   = afIn(i,p->t,0);
    float   *pfPlay = afPlayIn(i,p->t,0); 
    int32_t *pnRaw  = anIn(i,p->t,0);
    for (int n=0; n<p->Block; n++)
    {
        *pnRaw = *data;
        *pfIn  = (float)*data / 2147483648.0F;
        pfIn++; pnRaw++; data+=stride;
    }
}

void ThreadedDSP::Output(int o, int32_t *data, int stride) 
{ 
    if (!bOwner || o<0 || o>=p->O || !p || !p->bRunning) return;
    if (data == 0)
    {
        return;
    }

    if (!*abOut(o, p->t))                                          // Check for a Raw override - channel specific
    {
        int32_t *pN  = anTemp;
        ippsConvert_32f32s_Sfs(afOut(o,p->t,0), pN, p->Block, ippRndZero, -31);
        for (int n=0; n<p->Block; n++) 
        {
            *data = *pN++;
            data += stride;
        }
    }
    else
    {
        memcpy(data, anOut(o,p->t,0), p->N*sizeof(int32_t));      // buffer
        memset(anOut(o,p->t,0), 0, p->N*sizeof(int32_t));
        *abOut(o, p->t) = 0;
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
    p->Load = (float)tmax / p->N * nRate;
    Chrono_Load().add((int)(p->Load*100+0.5));
}

void ThreadedDSP::Process(void)
{
    if (!p || !bOwner) return;
    Chrono_CallTime().time();

    if (p->Filt_Next != p->Filt_Active) 
    {
        if (Filt_Set[p->Filt_Active]==0) Filt_Set[p->Filt_Active] = new Filter[p->F]();
        if (Filt_Set[p->Filt_Next]==0)   Filt_Set[p->Filt_Next]   = new Filter[p->F]();
        for (int f=0; f<p->F; f++) { p->Filt[f].Update = false; Filt_Set[p->Filt_Next][f].Update = false; };    // Avoid stale or lagging updates
        memcpy(Filt_Set[p->Filt_Active], p->Filt, p->F*sizeof(Filter));
        memcpy(p->Filt, Filt_Set[p->Filt_Next], p->F*sizeof(Filter));
        p->Filt_Active = p->Filt_Next;
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
            ippsAdd_32f_I(afPlayIn(n,p->t,0), afIn(n,p->t,0),p->Block); // Add the PlayIn buffer
            ippsMaxAbs_32f(afIn(n,p->t,0),p->Block,&fMax);              // Get the peak
            memset(afPlayIn(n,p->t,0),0,p->Block*sizeof(float));        // Clear the PlayIn buffer
            ippsMove_32f(afIn(n,(p->t-1+p->M)%p->M,0),   afTemp,         p->Block);	// 
            ippsMove_32f(afIn(n, p->t             ,0),   afTemp+p->Block,p->Block);	// Move two blocks of audio into FFT buffer
            ippsFFTFwd_RToPerm_32f(afTemp, afX(n,p->t,0), pFFT, pFFTBuf);  // And do the FFT
            if (fMax>p->PeakIn[n]) p->PeakIn[n]=fMax; else p->PeakIn[n] = 0.99F*p->PeakIn[n];  // Bit of smoothing
        }
        lk.lock(); if (++done==nThreads) middle.notify_all(); else middle.wait(lk); lk.unlock();    // Resync

        if (!p->bRunning) break;
        for (int n=ID; n<p->O; n+=nThreads) 			                // Do this over all of the outputs
		{				
			ippsSet_32f (0.0F, afY(n,0), 2*p->N);  				        // Zero MAC buffer
			for (int f=0; f<p->F; f++) 
			{															//
				if (p->Filt[f].Out != n) continue;		        		// Only filters with output channel n
                UpdateFilter(f, afTemp, pFFT, pFFTBuf);                 // Update the filter - out wont change, but m may
                for (int m=0; m<p->Filt[f].BLen; m++) 					// Loop over the filter M blocks
				{
					int nB = (p->t - m + p->M)%p->M;       		        // and do the MAC
                    float r0=*afY(n,0), r1=*afY(n,1);                   // Keep this as a multiple of 4 for SIMD optimization
					ippsAddProduct_32fc((const Ipp32fc *)afX(p->Filt[f].In,nB,0),(const Ipp32fc *)afH(f,m,0), (Ipp32fc *)afY(n,0), p->N);
					*afY(n,0) = r0 + *afX(p->Filt[f].In,nB,0) * *afH(f,m,0);    // Replace the first bin (DC + Nyquist)
					*afY(n,1) = r1 + *afX(p->Filt[f].In,nB,1) * *afH(f,m,1);    // with the normal multiply
				}
			}
 //         ippsMove_32f(afX(n,p->t,0), afY(n,0), 2*p->N);                  // FFT LOOPBACK
            ippsFFTInv_PermToR_32f(afY(n,0), afTemp, pFFT,pFFTBuf);		    // Do the inverse FFT and stick it in the
			ippsMove_32f(afTemp+p->N, afy(n,0), p->N);                      // correct output
            ippsMove_32f  (afy(n,0),           afOut(n,p->t,0),p->Block);   // Get the latest filter output
            ippsAdd_32f_I (afPlayOut(n,p->t,0),afOut(n,p->t,0),p->Block);   // Add in the Play buffer output
//          ippsMove_32f(afIn(n,p->t,0),afOut(n,p->t,0),p->Block);          // INPUT AUDIO LOOPBACK
            ippsSet_32f(0, afPlayOut(n,p->t,0),                p->Block);   // Clear the Play buffer
            ippsMaxAbs_32f(afOut(n,p->t,0),p->Block,&fMax);                 // Get the peak
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

void ThreadedDSP::UpdateFilter(int f,  float* afTemp, IppsFFTSpec_R_32f *pFFT, Ipp8u *pFFTBuf)
{
    if (!p->Filt[f].Update) return;
    int length   = p->Filt[f].Length;

    int M = (length-1)/p->N + 1;
    if (M != p->Filt[f].BLen) 
    {
        if (p->Filt[f].H) free(p->Filt[f].H);
        p->Filt[f].H = (float *)calloc(p->N*M,2*sizeof(float));
    }

    float* pFilt = afT(f,0);
    int        m = 0;
	while (length > 0 && m < p->M)
	{
		float *pF = afTemp;
		int     n = p->N;
		ippsSet_32f(0, afTemp, 2*p->N);
		while (n > 0 && length > 0) { *pF++ = *pFilt++; length--; n--; };
 		ippsFFTFwd_RToPerm_32f(afTemp, afH(f,m,0), pFFT, pFFTBuf);
		m++;
    }
	p->Filt[f].BLen = m;
    p->Filt[f].Update = false;
}

void ThreadedDSP::GetFilter(int in, int out, int maxlen, float *data)
{
    if (!p || !p->bRunning) return;        
	if (in < 0 || in >= p->I || out < 0 || out >= p->O) return;	
    if (maxlen==0) return;
    int f, m=0;
    for (f = 0; f<p->F; f++) if (p->Filt[f].In == in && p->Filt[f].Out == out) break;          // Check if we have this in->out already
	if (f == p->F) for (f = 0; f<p->F; f++) if (p->Filt[f].Length == 0) break;                 // Find the first free slot
    if (f == p->F) return;
    memcpy(data, afT(f,0), std::min(maxlen,p->Filt[f].Length)*sizeof(float));
}

int  ThreadedDSP::GetFilterLength(int in, int out)
{
    if (!p || !p->bRunning) return 0;        
	if (in < 0 || in >= p->I || out < 0 || out >= p->O) return 0;	
    int f, m=0;
    for (f = 0; f<p->F; f++) if (p->Filt[f].In == in && p->Filt[f].Out == out) break;          // Check if we have this in->out already
	if (f == p->F) for (f = 0; f<p->F; f++) if (p->Filt[f].Length == 0) break;                 // Find the first free slot
    if (f == p->F) return 0;
    return p->Filt[f].Length;
}


bool ThreadedDSP::LoadFilter(int in, int out, int length, float *pFilt)
{
    if (!p || !p->bRunning) return false;        
	if (in < 0 || in >= p->I || out < 0 || out >= p->O) return false;	
    
    int f, m=0;
    for (f = 0; f<p->F; f++) 
        if (p->Filt[f].In == in && p->Filt[f].Out == out) break;        // Check if we have this in->out already
    if (f == p->F)                                                      // 
    {                                                                   // If we don't have it
        if (length == 0) return true;                                   //   And if it is being cleared this is OK                                
    	for (f = 0; f<p->F; f++)                                        //
            if (p->Filt[f].Length==0 && p->Filt[f].BLen==0) break;      //   Otherwise search for a new slot
        if (f == p->F) return false;                                    //   If none available then return false
    }                                                                   // At this point we have an existing or new slot
    if (length>p->MBlock) length = p->MBlock;
    if (length>0) memcpy(afT(f,0), pFilt, length*sizeof(float));
    assert(p->Filt[f].BLen == 0 || (p->Filt[f].In==in && p->Filt[f].Out==out));                 // Just to be sure - not swapping live in/out
    p->Filt[f].In     = in;                 
    p->Filt[f].Out    = out;
    p->Filt[f].Length = length;
    p->Filt[f].Update = true;
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

void ThreadedDSP::GetIn  (int c, int s, float* data, uint64_t T, int stride)
{
    if (!p || !p->bRunning || c<0 || c>p->I || s<0 || s>p->MBlock || !data || T<0) return;
    if (T==0) T = p->t;
    int pos = ((T%p->M)*p->Block - s + p->MBlock)%(p->MBlock);  // Back in time - note p->t = p->T mod p->M
    int chunk = p->MBlock - pos;                    // Amount to read before buffer end
    if (chunk>s) chunk = s;                         // Only need at most s
    memcpy(data, afIn(c,0,0)+pos,chunk*sizeof(float)); data+=chunk; s-=chunk;
    if (s==0) return;
    memcpy(data, afIn(c,0,0),    s    *sizeof(float));
}

void ThreadedDSP::GetOut (int c, int s, float *data, uint64_t T, int stride)
{
    if (!p || !p->bRunning || c<0 || c>p->O || s<0 || s>p->MBlock || !data || T<0) return;
    if (T==0) T = p->t;
    int pos = ((T%p->M)*p->Block - s + p->MBlock)%(p->MBlock);  // Back in time - note p->t = p->T mod p->M
    int chunk = p->MBlock - pos;                    // Amount to read before buffer end
    if (chunk>s) chunk = s;                         // Only need at most s
    memcpy(data, afOut(c,0,0)+pos,chunk*sizeof(float)); data+=chunk; s-=chunk;
    if (s==0) return;
    memcpy(data, afOut(c,0,0),    s    *sizeof(float));
}

void ThreadedDSP::PlayIn (int c, int s, float* data, uint64_t T, int stride)
{
    if (!p || !p->bRunning || c<0 || c>p->I || s<0 || s>p->MBlock || !data || T<0) return;
    if (T==0) T = p->t;
    int chunk = p->Block*(p->M-(T%p->M));           // Current spot is our T modulus M
    if (chunk>s) chunk = s;                         // Only need to write up to s
    memcpy(afPlayIn(c,T%p->M,0), data, chunk*sizeof(float)); data+=chunk; s-=chunk;
    if (s==0) return;
    memcpy(afPlayIn(c,0    ,0), data, s    *sizeof(float));
}

void ThreadedDSP::PlayOut(int c, int s, float* data, uint64_t T, int stride)
{
    if (!p || !p->bRunning || c<0 || c>p->O || s<0 || s>p->MBlock || !data || T<0) return;
    if (T==0) T = p->t;
    int chunk = p->Block*(p->M-(T%p->M));           // Current spot is our T modulus M
    if (chunk>s) chunk = s;                         // Only need to write up to s
    memcpy(afPlayOut(c,T%p->M,0), data, chunk*sizeof(float)); data+=chunk; s-=chunk;
    if (s==0) return;
    memcpy(afPlayOut(c,0    ,0), data, s    *sizeof(float));
}

void ThreadedDSP::GetRaw  (int c, int s, int32_t* data, uint64_t T, int stride)
{
    if (!p || !p->bRunning || c<0 || c>p->I || s<0 || s>p->MBlock || !data || T<0) return;
    if (T==0) T = p->t;
    int pos = ((T%p->M)*p->Block - s + p->MBlock)%(p->MBlock);  // Back in time - note p->t = p->T mod p->M
    int chunk = p->MBlock - pos;                    // Amount to read before buffer end
    if (chunk>s) chunk = s;                         // Only need at most s
    memcpy(data, anIn(c,0,0)+pos,chunk*sizeof(int32_t)); data+=chunk; s-=chunk;
    if (s==0) return;
    memcpy(data, anIn(c,0,0),    s    *sizeof(int32_t));
}

void ThreadedDSP::PlayRaw(int c, int s, int32_t* data, uint64_t T, int stride)
{
    if (!p || !p->bRunning || c<0 || c>p->O || s<0 || s>p->MBlock || !data || T<0) return;
    if (T==0) T = p->t;
    int chunk = p->Block*(p->M-(T%p->M));           // Current spot is our T modulus M
    if (chunk>s) chunk = s;                         // Only need to write up to s
    for (int m=0; m<=(s-1)/p->N; m++) *abOut(c, (T+m)%p->M) = 1;        // Set the flag for output overwrite
    memcpy(anOut(c,T%p->M,0), data, chunk*sizeof(float)); data+=chunk; s-=chunk;
    if (s==0) return;
    memcpy(anOut(c,0    ,0), data, s    *sizeof(float));
}

