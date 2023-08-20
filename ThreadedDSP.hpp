#include <cstdint>
#include <thread>
#include <mutex>
#include <atomic>
#include <shared_mutex>
#include <condition_variable>
#include <cstring>
#include <math.h>
#include <histogram.hpp>

using namespace DA;

#ifdef _WIN32 
#include "../extern/ipp/include/ipp.h"
#else
#ifndef __arm__
//#include "../extern/ipp/include/ipp.h"
//#else 
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

#endif
#endif



#pragma warning( push )
#pragma warning( disable : 26451 )

class ThreadedDSP
{
public:
    static const int MaxChannels = 256;
    static const int MaxBlock    = 1024;
    static const int MaxLength   = 262144;
    static const int MaxThreads  = 32;
    static const int MaxFilters  = 8192;
    static const int MaxChronos  = MaxThreads+2+8;    // Spare chronos that can be shared

private: 
    struct Filter
    {
        int32_t     In;             // Input channel for this FIR
        int32_t     Out;            // Output channel
        int32_t     Length;         // Length of the filter in time domain taps
        int32_t     BLen;           // Length in integral blocks1 .. Blocks for the filter
        char        Name[64];       // A name for information or lookup
        int32_t     Update;         // Flag set to resynthesize from h -> H
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
        float       PeakIn[MaxChannels];    // Input peak values
        float       PeakOut[MaxChannels];   // Output peak values
        Filter      Filt[MaxFilters];       // Filter details
        float       Load;                   // An estimate of the load - DSP time / Call time


        histogram_t Chronos[MaxChronos];    // Shared chronos for CallTime, Load, ThreadExecTime and some spare
        int32_t     Threads;                //
        int32_t     Latency;                // Not used in this directly, but IO engine may set it
        float       Data[1];                // Member used to calculate the size of the Map structure prior to the audio data

/*      float       afIn[I][M][N];          // The Input time domain data           I*M*N       Data
        float       afOut[O][M][N];         // The Output time domain data          O*M*N       Data + ( I     )*M*N
        float       afPlayIn[I][M][N];      // Audio to play into the inputs        I*M*N       Data + ( I +  O)*M*N
        float       afPlayOut[O][M][N];     // Audio to play into the outputs       O*M*N       Data + (2I +  O)*M*N
        float       afX[I][M][2N];          // The buffer of F domain data          I*M*N*2     Data + (2I + 2O)*M*N
        float       afH[F][M][2N];          // The computed filter coefficients     F*M*N*2     Data + (4I + 2O)*M*N
        float       afT[F][M][N];           // The time domain filter coefficients  F*M*N       Data + (4I + 2O + 2F)*M*N
        float       afY[O][2N];             // The buffers for MAC for filter out   O*N*2       Data + (4I + 2O + 3F)*M*N
        float       afy[O][2N];             // Working buffer for output            O*N*2       Data + (4I + 2O + 3F)*M*N + 2*O*N
        int32_t     anIn[I][M][N]           // Bit perfect input Dante only         I*M*N       Data + (4I + 2O + 3F)*M*N + 4*O*N
        int32_t     anOut[O][M][N]          // Bit perfect output Dante override    O*M*N       Data + (5I + 2O + 3F)*M*N + 4*O*N
        int32_t     abOut[O][M]             // Flag set to raw ouput frame/channel  O*M         Data + (5I + 3O + 3F)*M*N + 4*O*N
        sizeof(Map) + ((5I + 3O + 3F)*M*N + 4*O*N + O*M)*sizeof(float)
*/
    };

    bool    bOwner;
    void*   hMem;
    char    sMem[256];
    size_t  nSize;
    Map*    p;         
    uint64_t LastT;
    int32_t FFTorder;               // Order of the FFT

public:
    float*  afIn    (int i, int m, int n) { return            p->Data +                            i   *p->MN + m*p->N   + n; };
    float*  afOut   (int o, int m, int n) { return            p->Data + (  p->I                  + o  )*p->MN + m*p->N   + n; };
    float*  afPlayIn(int i, int m, int n) { return            p->Data + (  p->I +   p->O         + i  )*p->MN + m*p->N   + n; };
    float*  afPlayOut(int o, int m, int n){ return            p->Data + (2*p->I +   p->O         + o  )*p->MN + m*p->N   + n; };
    float*  afX     (int i, int m, int n) { return            p->Data + (2*p->I + 2*p->O         + i*2)*p->MN + m*p->N*2 + n; };
    float*  afH     (int f, int m, int n) { return            p->Data + (4*p->I + 2*p->O         + f*2)*p->MN + m*p->N*2 + n; };
    float*  afT     (int f, int n)        { return            p->Data + (4*p->I + 2*p->O + 2*p->F+ f  )*p->MN            + n; };
    float*  afY     (int o, int n)        { return            p->Data + (4*p->I + 2*p->O + 3*p->F     )*p->MN +               + o*p->N*2 + n; };
    float*  afy     (int o, int n)        { return            p->Data + (4*p->I + 2*p->O + 3*p->F     )*p->MN + 2*p->O*p->N   + o*p->N*2 + n; };
    int32_t* anIn   (int i, int m, int n) { return (int32_t *)p->Data + (4*p->I + 2*p->O + 3*p->F+ i  )*p->MN + 4*p->O*p->N   + m*p->N   + n; };
    int32_t* anOut  (int o, int m, int n) { return (int32_t *)p->Data + (5*p->I + 2*p->O + 3*p->F+ o  )*p->MN + 4*p->O*p->N   + m*p->N   + n; };
    int32_t* abOut  (int o, int m       ) { return (int32_t *)p->Data + (5*p->I + 3*p->O + 3*p->F     )*p->MN + 4*p->O*p->N   + o*p->M   + m; };
    
    size_t   Size   (int I, int O, int M, int N, int F) { return sizeof(Map) +((5*I + 3*O + 3*F)*M*N + 4*O*N + O*M)*sizeof(float); };

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
    void GetIn  (int i, int s, float   *data, uint64_t T=0, int stride=1);     // Get the last s input samples of channel i 
    void GetOut (int o, int s, float   *data, uint64_t T=0, int stride=1);     // Get the last s output samples of channel o
    void PlayIn (int i, int s, float   *data, uint64_t T=0, int stride=1);     // Put s samples on the input queue for channel i
    void PlayOut(int o, int s, float   *data, uint64_t T=0, int stride=1);     // Push s samples on the output queue for channel o
    void GetRaw (int i, int s, int32_t *data, uint64_t T=0, int stride=1);     // Get the last s raw input sample of channel i
    void PlayRaw(int o, int s, int32_t *data, uint64_t T=0, int stride=1);     // Overwrite the next s output samples of channel o
    

    int     Inputs()    { if (!p) return 0; return p->I; };
    int     Outputs()   { if (!p) return 0; return p->O; };
    int     BlockSize() { if (!p) return 0; return p->Block; };
    int     Blocks()    { if (!p) return 0; return p->M; };
    int     BlockAt()   { if (!p) return 0; return p->t; };
    int     Filters()   { if (!p) return 0; int f=p->F; for (int n=0; n<p->F; n++) if (p->Filt[n].BLen) f--; return f; };
    uint64_t Count()     { if (!p) return 0; return *(volatile uint64_t *)&p->T; };
    bool    Running()   { if (!p) return 0; return *(volatile bool *)&p->bRunning; };
    bool    Stop()      { if (!p) return 0; p->bRunning=false; return true; };
    float   Load()      { if (!p) return 0; return p->Load; };
    bool    Owner()     { return bOwner; };
    int     Taps()      { if (!p) return 0; int nTaps=0; int f=p->F; for (int n=0; n<p->F; n++) nTaps+=p->Filt[n].BLen; return nTaps*p->Block; };
    int     Latency()   { if (!p) return 0; return p->Latency; };

    int     Threads()           { if (!p) return 0; return nThreads; };
    float   PeakIn(int n)       { if (!p || n<0 || n>p->I) return 0; return p->PeakIn[n]; };
    float   PeakOut(int n)      { if (!p || n<0 || n>p->O) return 0; return p->PeakOut[n]; };
    bool    LoadFilter      (int in, int out, int length=0, float *pFilt=0);
    void    GetFilter       (int in, int out, int maxlen, float *data);
    int     GetFilterLength (int in, int out);


    histogram_t *Chrono_CallTime(void)     { if (!p)                  return nullptr; return &p->Chronos[0]; };
    histogram_t *Chrono_ExecTime(int n)    { if (!p || n>=MaxThreads) return nullptr; return &p->Chronos[2+n]; };
    histogram_t *Chrono_Load(void)         { if (!p)                  return nullptr; return &p->Chronos[1]; };
    histogram_t *Chrono_N(int n)           { if (!p || n>=MaxChronos-MaxThreads-2) return nullptr; return &p->Chronos[MaxThreads+2+n]; };
    void Chrono_Reset(void)        
    { if (!p) return; for (int n=0; n<MaxChronos; n++) Histogram(&p->Chronos[n]).reset(); }; 

private:
    const int nRate = 48000;            // Maybe add an API to change later
    float   *afTemp;                    // Local buffer for filter upload
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
    void UpdateFilter(int f, float* afTemp, IppsFFTSpec_R_32f *pFFT, Ipp8u *pFFTBuf);
    std::chrono::high_resolution_clock::time_point Start;
    void ThreadProc(int ID);
    std::thread     *pThreads[MaxThreads];
    int              nThreads;
};