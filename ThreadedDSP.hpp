#include <cstdint>
#include <thread>
#include <mutex>
#include <atomic>
#include <shared_mutex>
#include <condition_variable>
#include <cstring>
#include <math.h>
#include <histogram.hpp>

using namespace DAES67;

#ifdef _WIN32 
#include "../extern/ipp/include/ipp.h"
#else
#if defined(__arm__) || defined(__aarch64__)
#include <arm_neon.h>
#include "../Ne10/inc/NE10.h"
#define ippStsNoErr          0
#define IPP_FFT_DIV_INV_BY_N 2
#define ippAlgHintNone       0
#define ippRndZero           0

typedef float    Ipp32f, Ipp32fc;
typedef int32_t  Ipp32s;
typedef uint8_t  Ipp8u;
typedef int      IppStatus;
typedef void     IppsFFTSpec_R_32f;
typedef int      IppRoundMode;
typedef int      IppHintAlgorithm;
inline Ipp8u*    ippsMalloc_8u         (int len) { if (len==0) return 0; else return (Ipp8u*)calloc(len,1); };
inline void      ippsFree              (void *p) { if (p) free(p); };
inline IppStatus ippsSet_32f           (Ipp32f val, Ipp32f* pDst, int len)                                                      { memset (pDst, 0,    len*sizeof(Ipp32f)); return ippStsNoErr; };
inline IppStatus ippsMove_32f          (const Ipp32f* pSrc, Ipp32f* pDst, int len)                                              { memmove(pDst, pSrc, len*sizeof(Ipp32f)); return ippStsNoErr; };
inline IppStatus ippsAdd_32f_I         (const Ipp32f* pSrc, Ipp32f* pSrcDst, int len)                                           { while(len) { vst1q_f32(pSrcDst,vaddq_f32(vld1q_f32(pSrcDst),vld1q_f32(pSrc))); pSrcDst+=4; pSrc+=4; len-=4; }; return ippStsNoErr; }
inline IppStatus ippsMaxAbs_32f        (const Ipp32f* pSrc, int len, Ipp32f* pMaxAbs)                                           { *pMaxAbs=0; while(len--) { if (fabs(*pSrc)>*pMaxAbs) *pMaxAbs=*pSrc; pSrc++; }; return ippStsNoErr; };
inline IppStatus ippsAddProduct_32fc   (const Ipp32fc* pSrc1, const Ipp32fc* pSrc2, Ipp32fc* pSrcDst, int len)                  
{   
    assert(len%4==0);
    float32x4x2_t A,B,C;
    while (len>=4)
    {
        A = vld2q_f32(pSrcDst);    B = vld2q_f32(pSrc1);    C = vld2q_f32(pSrc2);
        A.val[0] = vmlaq_f32(A.val[0],B.val[0],C.val[0]);
        A.val[0] = vmlsq_f32(A.val[0],B.val[1],C.val[1]);
        A.val[1] = vmlaq_f32(A.val[1],B.val[0],C.val[1]);
        A.val[1] = vmlaq_f32(A.val[1],B.val[1],C.val[0]);
        vst2q_f32(pSrcDst,A);
        pSrcDst+=8; pSrc1+=8; pSrc2+=8; len-=4;
    }
    return ippStsNoErr;
}
inline IppStatus ippsFFTGetSize_R_32f  (int order, int flag, IppHintAlgorithm hint, int*pSpecSize, int* pSpecBufferSize, int* pBufferSize)           { *pSpecSize=0; pSpecBufferSize=0; *pBufferSize = (int)sizeof(float32_t)*((1<<order)+2); return ippStsNoErr; };
inline IppStatus ippsFFTInit_R_32f     (IppsFFTSpec_R_32f** ppFFTSpec, int order, int flag, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer) { *ppFFTSpec = (IppsFFTSpec_R_32f *)ne10_fft_alloc_r2c_float32(1<<order); return ippStsNoErr; };
inline IppStatus ippSetDenormAreZeros  (int value)                                                                                                   { return ippStsNoErr; };
inline IppStatus ippsFFTFwd_RToPerm_32f(Ipp32f* pSrc, Ipp32f* pDst, const IppsFFTSpec_R_32f* pFFTSpec, Ipp8u* pBuffer)    
{ 
    ne10_fft_r2c_cfg_float32_t FFT = (ne10_fft_r2c_cfg_float32_t)pFFTSpec;
    float32_t *Temp = (float32_t *)pBuffer;
    int32_t   N = FFT->nfft;
    ne10_fft_r2c_1d_float32_c((ne10_fft_cpx_float32_t *)Temp, pSrc, FFT);
    pDst[0]=Temp[0]; pDst[1]=Temp[N]; memmove(pDst+2,Temp+2,(N-2)*sizeof(Ipp32f));
    return ippStsNoErr; 
};

inline IppStatus ippsFFTInv_PermToR_32f(const Ipp32f* pSrc, Ipp32f* pDst, const IppsFFTSpec_R_32f* pFFTSpec, Ipp8u* pBuffer)
{
    ne10_fft_r2c_cfg_float32_t FFT = (ne10_fft_r2c_cfg_float32_t)pFFTSpec;
    float32_t *Temp = (float32_t *)pBuffer;
    int32_t   N = FFT->nfft;
    Temp[0]=pSrc[0]; Temp[1]=0; Temp[N]=pSrc[1]; Temp[N+1]=0; memmove(Temp+2,(void *)(pSrc+2),(N-2)*sizeof(Ipp32f));
    ne10_fft_c2r_1d_float32_c(pDst, (ne10_fft_cpx_float32_t *)Temp, FFT);
    return ippStsNoErr; 
};

inline IppStatus ippsConvert_32f32s_Sfs(const Ipp32f* pSrc, Ipp32s* pDst, int len, IppRoundMode rndMode, int scaleFactor)                            
{ 
    const float fClamp = (float)0x7FFF0000L;
    float       fScale = (float)(1L<<(-scaleFactor));
    for (int n=0; n<len; n++) 
    { 
        float fTmp = fScale * *pSrc++; 
        if      (fTmp >  fClamp) fTmp =  fClamp;
        else if (fTmp < -fClamp) fTmp = -fClamp;;
        *pDst++=(Ipp32s)(fTmp + 0.5F);
    }
    return ippStsNoErr; 
};
//  COULD OPTIMIZE THIS { while(len) { vst1q_f32(pSrcDst,vaddq_f32(vld1q_f32(pSrcDst),vld1q_f32(pSrc))); pSrcDst+=4; pSrc+=4; len-=4; }; return ippStsNoErr; }

#else 
#include "ipp.h"
#endif
#endif



#pragma warning( push )
#pragma warning( disable : 26451 )

class ThreadedDSP
{
public:
    static const int MaxChannels = 512;
    static const int MaxBlock    = 1024;
    static const int MaxLength   = 262144;
    static const int MaxThreads  = 32;
    static const int MaxFilters  = 512;
    static const int MaxChronos  = MaxThreads+2+8;    // Spare chronos that can be shared
    static const int MaxSets     = 1512;
    static const int MaxGroups   = 2;

private: 
    struct Filter
    {
        int32_t     In;             // Input channel for this FIR
        int32_t     Out;            // Output channel
        int32_t     Length;         // Length of the filter in time domain taps
        int32_t     BLen;           // Length in integral blocks1 .. Blocks for the filter
        int32_t     Update;         // Flag set to resynthesize from h -> H
        float*      H = nullptr;    // The frequency domain coefficients (only valid poitner in processing thread)
    };

    struct Map
    {
        bool        bRunning;               // Set by owner to false when exiting
        int32_t     Block;                  // Number of samples in each call
        int32_t     N;                      // Block + XFade/2
        int32_t     M;                      // Number of Blocks in filter and IO buffer
        int32_t     I;                      // Number of input channels
        int32_t     O;                      // Number of output channels
        int32_t     F;                      // Numbwe of filters
        int32_t     t;                      // Current block in the M cycle
        int64_t     T;                      // Number of blocks processed
        int32_t     MBlock;
        int32_t     MN;
        float       GainIn[MaxChannels];    // Input gains
        float       GainOut[MaxChannels];   // Output gains
        float       PeakIn[MaxChannels];    // Input peak values
        float       PeakOut[MaxChannels];   // Output peak values
        Filter      Filt[MaxGroups][MaxFilters];       // Current active Filter details
        int32_t     Filt_Active[MaxGroups];            // Active filters
        int32_t     Filt_Next[MaxGroups];              // Next filters to swap in asap
        bool        ClearAll;               // Flag to clear all filters
        float       Load;                   // An estimate of the load - DSP time / Call time


        Histogram   Chronos[MaxChronos];    // Shared chronos for CallTime, Load, ThreadExecTime and some spare
        int32_t     Threads;                //
        int32_t     Latency;                // Not used in this directly, but IO engine may set it
        float       Data[1];                // Member used to calculate the size of the Map structure prior to the audio data

/*      float       afIn[I][M][N];          // The Input time domain data           I*2*N       Data
        float       afOut[O][M][N];         // The Output time domain data          O*2*N       Data +  I     *2*N
        float       afT[F][M][N];           // The time domain filter coefficients  F*MBlock    Data + (I + O)*2*N  
        sizeof(Map) + ((I + O)*2*N F*MBlock)*sizeof(float)
*/
    };

    bool    bOwner;
    void*   hMem;
    char    sMem[256];
    size_t  nSize;
    Map*    p;         
    uint64_t LastT;

private:   
    float* pfX;            // The buffer of F domain data        I*M*N*2
    float* pfY;            // The buffers for MAC for filter out   O*N*2
    float* pW;             // The window for the output overlap    2*N - Block

public:
    float*  afIn    (int i, int n)        { return p->Data +      i         *2*p->N   + n; };
    float*  afOut   (int o, int n)        { return p->Data + ( p->I +    o )*2*p->N   + n; };
    float*  afT     (int f, int n)        { return p->Data + ( p->I + p->O )*2*p->N   + f*p->MBlock + n; };
    size_t   Size   (int I, int O, int M, int N, int F, int B) { return sizeof(Map) +((I + O)*2*N + F*M*B)*sizeof(float); };


    float*  afX     (int i, int m, int n) { return pfX + i*2*p->MN + m*p->N*2 + n; };
    float*  afH     (int f, int m, int n, int g=0) { static float zero[2*MaxBlock];  if (p->Filt[g][f].H==nullptr) return zero; else return p->Filt[g][f].H + m*p->N*2 + n; };
    float*  afY     (int o, int n)        { return pfY + o*p->N*2 + n; };
    

public:
    ThreadedDSP();          
   ~ThreadedDSP();

    bool Create(int Block, int Blocks,                  // Create new Threaded DSP with block size and number of blocks (max filter length)
                int In,      int Out,                   // Set the number of input and output channels
                int Filters,                            // Set the maximum number of filters
                int XFade = 0,                          // Set the size of the cross fade - FFT blocks are 2*Block + XFade
                const char* sShared = "ThreadedDSP",    // The string name for the shared memory
                int nThreads = 4,                       // The number of threads
                int nLatency = 0 );                     // An optional latency value in samples by creator engine
    bool Attach(const char* sShared = "ThreadedDSP");   // Attach this class to an existing shared memory                                                                
    void Destroy();

    void Input (int n, int32_t* data, int stride=1);
    void Input (int n, float*   data, int stride=1);
    void Output(int n, int32_t* data, int stride=1);
    void Output(int n, float*   data, int stride=1);
    const char* SHM_Name() { return sMem; };

    void Finish(void);
    void Process();

    // If T=0 it will return the data offset from the current p->t
    // Otherwise it will read/write as if p->T==T, which is eqivalent to p->t = T%p->M
    // In the case of T!=0, reading or writing a full set of Block*M samples will cause a wrap.
    // These are non blocking, so if using T==0, use Spin to ensure data consistency
    uint64_t Spin   (uint64_t T=0, int timeoutms=500, bool spin=false);           // Spin a hard or sleep cycle until just after last block
    uint64_t SpinT  () { return LastT; };                                    // The frame at exit of last Spin (same as return from last Spin)

    int     Inputs()    { if (!p) return 0; return p->I; };
    int     Outputs()   { if (!p) return 0; return p->O; };
    int     BlockSize() { if (!p) return 0; return p->Block; };
    int     FFTSize()   { if (!p) return 0; return p->N; };
    int     CrossFade() { if (!p) return 0; return 2*(p->N - p->Block); };
    int     Blocks()    { if (!p) return 0; return p->M; };
    int     BlockAt()   { if (!p) return 0; return p->t; };
    uint64_t Count()     { if (!p) return 0; return *(volatile uint64_t *)&p->T; };
    bool    Running()   { if (!p) return 0; return *(volatile bool *)&p->bRunning; };
    bool    Stop()      { if (!p) return 0; p->bRunning=false; return true; };
    float   Load()      { if (!p) return 0; return p->Load; };
    bool    Owner()     { return bOwner; };
    int     Filters(int g=0)   { if (!p || g<0 || g>=MaxGroups) return 0; int f=0;     for (int n=0; n<p->F; n++) if (p->Filt[g][n].BLen) f++; return f; };
    int     Taps(int g=0)      { if (!p || g<0 || g>=MaxGroups) return 0; int nTaps=0; for (int n=0; n<p->F; n++) nTaps+=p->Filt[g][n].BLen; return nTaps*p->Block; };
    int     Latency()   { if (!p) return 0; return p->Latency; };

    int     Threads()           { if (!p) return 0; return nThreads; };
    float   PeakIn(int n)       { if (!p || n<0 || n>p->I) return 0; return p->PeakIn[n]; };
    float   PeakOut(int n)      { if (!p || n<0 || n>p->O) return 0; return p->PeakOut[n]; };
    bool    LoadFilter      (int in, int out, int length=0, float *pFilt=0, int group=0);

    float   GetGainIn(int n)  { if (!p || n<0 || n>p->I) return 0; return p->GainIn[n]; };
    float   GetGainOut(int n) { if (!p || n<0 || n>p->O) return 0; return p->GainOut[n]; };
    void    SetGainIn(int n, float f)  { if (!p || n<0 || n>p->I) return; p->GainIn[n] = f;  };
    void    SetGainOut(int n, float f) { if (!p || n<0 || n>p->O) return; p->GainOut[n] = f; };
    void    SetGainsIn(float* fIn)     { if (!p) return; for (int n=0; n<p->I; n++) p->GainIn[n] = fIn[n]; };
    void    SetGainsOut(float* fOut)   { if (!p) return; for (int n=0; n<p->O; n++) p->GainOut[n] = fOut[n]; };

    int     GetFilterSet    (int g=0)      { if (g<0 || g>MaxGroups) return 0; return p->Filt_Active[g]; };
    int     SetFilterSet    (int n, int g=0) 
    { 
        if (p==0 || n<0 || n>=MaxSets) return 0;
        if (g<0 || g>=MaxGroups) return 0;
        if (p->Filt_Active[g] == n) return n;
        int old = p->Filt_Active[g]; p->Filt_Next[g] = n; 
        return old; 
    };
    void   ClearFilters() { if (p==0 ) return; p->ClearAll = true; };

    Histogram& Chrono_CallTime(void)     { return p->Chronos[0]; };
    Histogram& Chrono_ExecTime(int n)    { return p->Chronos[2+n]; };
    Histogram& Chrono_Load(void)         { return p->Chronos[1]; };
    Histogram& Chrono_N(int n)           { if (n<MaxChronos-MaxThreads-2) return p->Chronos[MaxThreads+2+n]; else return p->Chronos[MaxThreads+2]; };
    void Chrono_Reset(void)        
    { if (!p) return; for (int n=0; n<MaxChronos; n++) p->Chronos[n].reset(); }; 

private:
    const int nRate = 48000;            // Maybe add an API to change later
    int32_t *anTemp;                    // Local buffer for data download
	Ipp8u *pFFTSpec, *pFFTInit, *pFFTBuf;
	IppsFFTSpec_R_32f *pFFT;
    std::mutex              mtx;        // For scheduling threads in the owner
    std::condition_variable start;      // Condition used to kick off all threads
    std::condition_variable middle;     // Condition used to signal that all threads are done
    std::condition_variable finished;   // Condition used to signal that all threads are done
    std::atomic<int>        done;       // Count of threads as they complete each stage
    std::atomic<bool>       processing; // Set between Process and Finish
    std::atomic<bool>       loading;    // Set when loading filters
    void UpdateFilter(int f, float* afTemp, IppsFFTSpec_R_32f *pFFT, Ipp8u *pFFTBuf, int group=0);
    std::chrono::high_resolution_clock::time_point Start;
    void ThreadProc(int ID);
    std::thread     *pThreads[MaxThreads];
    int              nThreads;

    Filter *Filt_Set[MaxGroups][MaxSets] = {};
};