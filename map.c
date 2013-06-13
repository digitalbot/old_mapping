/*
 * 到達時間差を用いた近接4点法の推定可能点分布を作成する
 * usage: gcc map.c -o mapping_calc; ./mapping_calc
 */

/*
 *
 * [0]   [-15] = -0.000156
 * [1]   [-14] = -0.000146
 * [2]   [-13] = -0.000135
 * [3]   [-12] = -0.000125
 * [4]   [-11] = -0.000115
 * [5]   [-10] = -0.000104
 * [6]   [-9]  = -0.000094
 * [7]   [-8]  = -0.000083
 * [8]   [-7]  = -0.000073
 * [9]   [-6]  = -0.000063
 * [10]  [-5]  = -0.000052
 * [11]  [-4]  = -0.000042
 * [12]  [-3]  = -0.000031
 * [13]  [-2]  = -0.000021
 * [14]  [-1]  = -0.000010
 * [15]  [0]   =  0.000000
 * [16]  [1]   =  0.000010
 * [17]  [2]   =  0.000021
 * [18]  [3]   =  0.000031
 * [19]  [4]   =  0.000042
 * [20]  [5]   =  0.000052
 * [21]  [6]   =  0.000063
 * [22]  [7]   =  0.000073
 * [23]  [8]   =  0.000083
 * [24]  [9]   =  0.000094
 * [25]  [10]  =  0.000104
 * [26]  [11]  =  0.000115
 * [27]  [12]  =  0.000125
 * [28]  [13]  =  0.000135
 * [29]  [14]  =  0.000146
 * [30]  [15]  =  0.000156
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

const double DIST         = 0.05;
const double SONIC        = 340.0;
const double MAX_LAG      = 0.0001470588235;

double FS = 96000.0;           // 標本化周波数
double TS = 0.0000104166667;   // 標本化周期

typedef struct {
    double x;
    double y;
    double z;
} st_answer;

st_answer calc_answer(double a_lag, double b_lag, double c_lag);
void pre_calc(double *list, int limit, int  range, double mag);


int main(int argc, char *argv[]) {

    char ratio[5];
    printf("type number [1~32]. (32 is very very heavy...)\n");
    fgets(ratio, 5, stdin);
    double mag = atof(ratio);

    int limit_sample    = (MAX_LAG / (TS / mag)) + 1.5;
    int range           = limit_sample * 2;

    printf("range is %d\n", range);

    /* these are dist lag */
    double *lag_list = calloc(range + 1, sizeof(double));

    int i;
    for (i=0; i<=range; i++) {
        lag_list[i] = (i - limit_sample) * (TS / mag);
        printf("%d %.10f\n", i - limit_sample, lag_list[i]);
    }

    pre_calc(lag_list, limit_sample, range, mag);

    printf("DONE\n");
    printf("ratio %f, range %d\n", mag, range);
    printf("max index: %d\n", (int)pow((range + 1), 3));
    free(lag_list);
    return 0;
}

void pre_calc(double *list, int limit, int range, double mag) {
    char f_name[100];
    sprintf(f_name, "mapping_%f.csv", mag);
    FILE *fp = fopen(f_name, "w");
    if (fp == NULL) {
        printf("ERROR!\n");
        exit(-1);
    }
    int j, k, l;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (j=0; j<=range; j++) {
        for (k=0; k<=range; k++) {
            for (l=0; l<=range; l++) {
                st_answer answers = calc_answer(list[j], list[k], list[l]);
                //printf("%d:%d:%d -> %d, %f, %f, %f\n",
                //       j, k, l, (j * (range + 1) * (range + 1)) +(k * (range + 1)) + l,
                //       answers.x, answers.y, answers.z);
                if (answers.x == 987654321.0) {
                }
                else{
                    fprintf(fp, "%d, %.13f, %.13f, %.13f\n",
                            (j * (range + 1) * (range + 1)) + (k * (range + 1)) + l,
                            answers.x, answers.y, answers.z);
                }
            }
            printf(".\n");
        }
        printf("_____\n");
    }
    fclose(fp);
    printf("Output was named %s\n", f_name);
    return;
}

st_answer calc_answer(double a_lag, double b_lag, double c_lag) {

    st_answer answers;
    answers.x = 0;
    answers.y = 0;
    answers.z = 0;

    double drA = SONIC * a_lag;
    double drB = SONIC * b_lag;
    double drC = SONIC * c_lag;
    /* powered */
    double p2_drA  = pow(drA, 2.0);
    double p2_drB  = pow(drB, 2.0);
    double p2_drC  = pow(drC, 2.0);
    double p2_dist = pow(DIST, 2.0);

    /* put BIG A, B, C to tmpOne, tmpTwo, tmpThree */
    double tmpOne   = 4.0 * (+ (3.0 * (+ p2_drA
                                       + p2_drB
                                       + p2_drC)
                                )
                             - 2.0 * (+ (drA * drB)
                                      + (drA * drC)
                                      + (drB * drC)
                                      + p2_dist)
                             );

    double tmpTwo   = 4.0 * (- 3.0 * (+ pow(drA, 3.0)
                                      + pow(drB, 3.0)
                                      + pow(drC, 3.0)
                                      )
                             + (p2_drA * (drB + drC))
                             + (p2_drB * (drA + drC))
                             + (p2_drC * (drA + drB))
                             + (p2_dist * (drA + drB + drC))
                             );

    double tmpThree = 1.0 * (+ (3.0 * (+ pow(DIST, 4.0)
                                       + pow(drA, 4.0)
                                       + pow(drB, 4.0)
                                       + pow(drC, 4.0)
                                       )
                                )
                             - ((2.0 * p2_dist) * (p2_drA + p2_drB + p2_drC))
                             - 2.0 * (+ (p2_drA * p2_drB)
                                      + (p2_drC * (p2_drA + p2_drB)
                                         )
                                      )
                             );
    double rO;
    if (tmpOne == 0.0) {
        rO = - 1.0 * (tmpThree / tmpTwo);
    }
    else {
        rO = - 1.0 * (tmpTwo / (2.0 * tmpOne))
             + (sqrt(pow(tmpTwo, 2.0)
                     - (4.0 * tmpOne * tmpThree)
                     )
                / fabs(2.0 * tmpOne));
    }
    if (isnan(rO)) {
        answers.x = 987654321.0;
        answers.y = 987654321.0;
        answers.z = 987654321.0;
        return answers;
    }
    double rA = rO - drA;
    double rB = rO - drB;
    double rC = rO - drC;
    double p2_rO = pow(rO, 2);
    double p2_rA = pow(rA, 2);
    double p2_rB = pow(rB, 2);
    double p2_rC = pow(rC, 2);

    double tmpX = (p2_rA - p2_dist - p2_rO) / (2 * DIST);
    double tmpY = (+ (2 * p2_rB)
                   - p2_dist
                   - p2_rO
                   - p2_rA
                   ) / (2 * sqrt(3) * DIST);
    double tmpZ = (+ (3 * p2_rC)
                   - p2_dist
                   - p2_rO
                   - p2_rA
                   - p2_rB
                   ) / (2 * sqrt(6) * DIST);

    answers.x = tmpX;
    answers.y = tmpY;
    answers.z = tmpZ;
    return answers;
}
