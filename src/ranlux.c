#include <stdio.h>
#include <stdlib.h>

#define MAXLEV 4
#define LXDFLT 3
#define TWOP12 4096
#define IGIGA 1000000000
#define JSDFLT 314159265
#define ITWO24 (1 << 24)
#define ICONS 2147483563

// 定义结构体来保存状态
typedef struct {
    int LUXLEV;
    int NOTYET;
    int IN24;
    int KOUNT;
    int MKOUNT;
    int I24;
    int J24;
    double CARRY;
    double SEEDS[24];
    double TWOM24;
    double TWOM12;
    int NSKIP;
    int NDSKIP[MAXLEV + 1];
    int NEXT[24];
    int INSEED;
} RANLUX_STATE;

// 初始化全局变量
RANLUX_STATE state = {
   .LUXLEV = LXDFLT,
   .NOTYET = 1,
   .IN24 = 0,
   .KOUNT = 0,
   .MKOUNT = 0,
   .I24 = 24,
   .J24 = 10,
   .CARRY = 0,
   .TWOM24 = 1,
   .TWOM12 = TWOM24 * TWOP12,
   .NSKIP = state.NDSKIP[LXDFLT],
};

void RANLUX(double *RVEC, int LENV) {
    if (state.NOTYET) {
        state.NOTYET = 0;
        int JSEED = JSDFLT;
        state.INSEED = JSEED;
        printf(" RANLUX DEFAULT INITIALIZATION: %d\n", JSEED);
        state.LUXLEV = LXDFLT;
        state.NSKIP = state.NDSKIP[state.LUXLEV];
        int LP = state.NSKIP + 24;
        state.IN24 = 0;
        state.KOUNT = 0;
        state.MKOUNT = 0;
        printf(" RANLUX DEFAULT LUXURY LEVEL =  %d      p =%d\n", state.LUXLEV, LP);
        state.TWOM24 = 1;
        for (int i = 1; i <= 24; i++) {
            state.TWOM24 *= 0.5;
        }
        int K;
        for (int i = 1; i <= 24; i++) {
            K = JSEED / 53668;
            JSEED = 40014 * (JSEED - K * 53668) - K * 12211;
            if (JSEED < 0) JSEED += ICONS;
            int ISEEDS = JSEED % ITWO24;
            state.SEEDS[i - 1] = (double)ISEEDS * state.TWOM24;
        }
        for (int i = 1; i <= 24; i++) {
            state.NEXT[i - 1] = i - 1;
        }
        state.NEXT[0] = 24;
        state.I24 = 24;
        state.J24 = 10;
        state.CARRY = 0;
        if (state.SEEDS[23] == 0) state.CARRY = state.TWOM24;
    }

    for (int IVEC = 0; IVEC < LENV; IVEC++) {
        double UNI = state.SEEDS[state.J24] - state.SEEDS[state.I24] - state.CARRY;
        if (UNI < 0) {
            UNI += 1.0;
            state.CARRY = state.TWOM24;
        } else {
            state.CARRY = 0;
        }
        state.SEEDS[state.I24] = UNI;
        state.I24 = state.NEXT[state.I24];
        state.J24 = state.NEXT[state.J24];
        RVEC[IVEC] = UNI;
        if (UNI < state.TWOM12) {
            RVEC[IVEC] += state.TWOM24 * state.SEEDS[state.J24];
            if (RVEC[IVEC] == 0) RVEC[IVEC] = state.TWOM24 * state.TWOM24;
        }
        state.IN24++;
        if (state.IN24 == 24) {
            state.IN24 = 0;
            state.KOUNT += state.NSKIP;
            for (int ISK = 0; ISK < state.NSKIP; ISK++) {
                double UNI = state.SEEDS[state.J24] - state.SEEDS[state.I24] - state.CARRY;
                if (UNI < 0) {
                    UNI += 1.0;
                    state.CARRY = state.TWOM24;
                } else {
                    state.CARRY = 0;
                }
                state.SEEDS[state.I24] = UNI;
                state.I24 = state.NEXT[state.I24];
                state.J24 = state.NEXT[state.J24];
            }
        }
    }
    state.KOUNT += LENV;
    if (state.KOUNT >= IGIGA) {
        state.MKOUNT++;
        state.KOUNT -= IGIGA;
    }
}

void RLUXGO(int LUX, int INS, int K1, int K2) {
    if (LUX < 0) {
        state.LUXLEV = LXDFLT;
    } else if (LUX <= MAXLEV) {
        state.LUXLEV = LUX;
    } else if (LUX < 24 || LUX > 2000) {
        state.LUXLEV = MAXLEV;
        printf(" RANLUX ILLEGAL LUXURY RLUXGO: %d\n", LUX);
    } else {
        state.LUXLEV = LUX;
        for (int ILX = 0; ILX <= MAXLEV; ILX++) {
            if (LUX == state.NDSKIP[ILX] + 24) state.LUXLEV = ILX;
        }
    }
    if (state.LUXLEV <= MAXLEV) {
        state.NSKIP = state.NDSKIP[state.LUXLEV];
        printf(" RANLUX LUXURY LEVEL SET BY RLUXGO :%d     P=%d\n", state.LUXLEV, state.NSKIP + 24);
    } else {
        state.NSKIP = state.LUXLEV - 24;
        printf(" RANLUX P - VALUE SET BY RLUXGO TO:%d\n", state.LUXLEV);
    }
    state.IN24 = 0;
    if (INS < 0) {
        printf(" Illegal initialization by RLUXGO, negative input seed\n");
    }
    int JSEED;
    if (INS > 0) {
        JSEED = INS;
        printf(" RANLUX INITIALIZED BY RLUXGO FROM SEEDS %d %d %d\n", JSEED, K1, K2);
    } else {
        JSEED = JSDFLT;
        printf(" RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED\n");
    }
    state.INSEED = JSEED;
    state.NOTYET = 0;
    state.TWOM24 = 1;
    for (int i = 1; i <= 24; i++) {
        state.TWOM24 *= 0.5;
    }
    int K;
    for (int i = 1; i <= 24; i++) {
        K = JSEED / 53668;
        JSEED = 40014 * (JSEED - K * 53668) - K * 12211;
        if (JSEED < 0) JSEED += ICONS;
        int ISEEDS = JSEED % ITWO24;
        state.SEEDS[i - 1] = (double)ISEEDS * state.TWOM24;
    }
    for (int i = 1; i <= 24; i++) {
        state.NEXT[i - 1] = i - 1;
    }
    state.NEXT[0] = 24;
    state.I24 = 24;
    state.J24 = 10;
    state.CARRY = 0;
    if (state.SEEDS[23] == 0) state.CARRY = state.TWOM24;
    state.KOUNT = K1;
    state.MKOUNT = K2;
    if (K1 + K2!= 0) {
        for (int IOUTER = 1; IOUTER <= K2 + 1; IOUTER++) {
            int INNER = IGIGA;
            if (IOUTER == K2 + 1) INNER = K1;
            for (int ISK = 0; ISK < INNER; ISK++) {
                double UNI = state.SEEDS[state.J24] - state.SEEDS[state.I24] - state.CARRY;
                if (UNI < 0) {
                    UNI += 1.0;
                    state.CARRY = state.TWOM24;
                } else {
                    state.CARRY = 0;
                }
                state.SEEDS[state.I24] = UNI;
                state.I24 = state.NEXT[state.I24];
                state.J24 = state.NEXT[state.J24];
            }
        }
        state.IN24 = state.KOUNT % (state.NSKIP + 24);
        if (state.MKOUNT > 0) {
            int IZIP = IGIGA % (state.NSKIP + 24);
            int IZIP2 = state.MKOUNT * IZIP + state.IN24;
            state.IN24 = IZIP2 % (state.NSKIP + 24);
        }
        if (state.IN24 > 23) {
            printf("  Error in RESTARTING with RLUXGO:\n  The values %d %d %d cannot occur at luxury level %d\n", INS, K1, K2, state.LUXLEV);
            state.IN24 = 0;
        }
    }
}

int main() {
    double random_vector[10];
    RANLUX(random_vector, 10);
    for (int i = 0; i < 10; i++) {
        printf("%lf\n", random_vector[i]);
    }
    return 0;
}