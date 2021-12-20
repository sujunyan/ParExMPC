#include "pempc.h"
#include "mpQP_data.h"
#include "MPC_problem.h"

void setup() {
  // setup
  Serial.begin(9600);
  Serial.println("begin");
}

void loop() {
    static double x0[] = {10,0};
    static double u0[nu];
    static double yin[y_size];
    static double lam_in[lam_size];
    static double yout[y_size];
    static double lam_out[lam_size];
    static size_t max_iter = 3;
    static double tol = 1e-3;
    size_t time_;
    pempc_aAx('n',1,A,x0,nx,nx,x0);
    pempc_aAxpy('n',1,B,u0,nx,nu,x0);
    time_ = micros();
    pempc_get_control(x0,max_iter,tol,yout,lam_out,yout,lam_out,u0);
    time_ = micros() -time_;
    Serial.println(time_/1000.0f);
    Serial.println(u0[0]);
    delay(1000);
}
