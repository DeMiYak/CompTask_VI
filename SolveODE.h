//
// Created by Пользователь on 17.05.2023.
//

#ifndef COMPTASK_VI_SOLVEODE_H
#define COMPTASK_VI_SOLVEODE_H

#include <boost/numeric/ublas/matrix.hpp>
#include"RPN.h"

namespace ublas = boost::numeric::ublas;

class SolveODE{
public:
    SolveODE() = default;

    SolveODE(const double &xCauchy, const double &yCauchy){
        _xCauchy = xCauchy;
        _yCauchy = yCauchy;
    }

    SolveODE(const double &xCauchy, const double &yCauchy, const double &step, const int &stepNum, const Formula &DifferentialEquation){
        _xCauchy = xCauchy;
        _yCauchy = yCauchy;

        if(GoodInput(step, stepNum)){
            _step = step;
            _stepNum = stepNum;
        }

        _AdamsExtrapolationMatrix.resize(stepNum + 2, 7);
        _DifferentialEquation = DifferentialEquation;

        this->HardcodeTaylorCoefficient(6);
        this->BuildTaylorSolution();
        this->InsertTaylorToAdamsMatrix();
    }

    void HardcodeTaylorCoefficient(const int &coefficientNum);

    void BuildTaylorSolution();

    double EvaluateTaylorSolution(double x);

    void PrintTaylorValues();

    void PrintPreciseValues(Formula &referenceFormula);

    void ErrorTaylorValues(Formula &referenceFormula);

    void PrintValues();


    void InsertTaylorToAdamsMatrix();

    void BuildAdamsExtrapolationMatrix();

    void PrintAdamsMatrix();

    void PrintExtrapolationValues();

    double EvaluateDifferentialEquation(double x, double y);


    ublas::vector<double> Runge_KuttaCoefficient(const double &x, const double &y);

    void PrintRunge_Kutta();


    void EulerMethod0();

    void EulerMethodI();

    void EulerMethodII();


    double GetXCauchy(){return _xCauchy;}

    double GetYCauchy(){return _yCauchy;}

    double GetStep(){return _step;}

    int GetStepNum(){return _stepNum;}

    Formula GetTaylorSolution(){return _TaylorSolution;}

    ublas::vector<double> GetTaylorCoefficient(){return _TaylorCoefficient;}

private:

    double _xCauchy = 0;

    double _yCauchy = 1;

    double _step = 0.1;

    int _stepNum = 10;

    Formula _TaylorSolution = Formula("");

    Formula _DifferentialEquation = Formula("");

    ublas::vector<double> _TaylorCoefficient = ublas::vector<double>(6, 0);

    ublas::matrix<double> _AdamsExtrapolationMatrix = ublas::matrix<double>(5, 7);

    bool GoodInput(const double &step, const double &stepNum){
        if(step > 0 && stepNum > 1) return true;
        return false;
    }
};

#endif //COMPTASK_VI_SOLVEODE_H
