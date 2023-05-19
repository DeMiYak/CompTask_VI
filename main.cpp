#include <iostream>
#include<iomanip>
#include"RPN.h"
#include"SolveODE.h"
#include"TaylorTask.h"
#include"AdamsTask.h"
#include"Runge-KuttaTask.h"
#include"EulerTask.h"

int main() {
    _wsetlocale(LC_ALL, L"russian");
    wcout.precision(16);
    cout.precision(16);

    TaylorTask();


    /*
     * On launch references out of boundary values in ODE.PrintAdamsMatrix();
     * When debugged or some other function is utilized before it
     * Prints well-reserved values.
     * Check compiler in this regard.
     */
    AdamsTask();
//
//    AdamsTask();

//    Runge_KuttaTask();

//    EulerTask();
}
