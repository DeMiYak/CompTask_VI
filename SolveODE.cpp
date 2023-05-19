//
// Created by Пользователь on 17.05.2023.
//

#include"SolveODE.h"

int factorial(int number) {
    int factorial = 1;
    for (int i = 2; i < number; ++i) {
        factorial *= i;
    }
    return factorial;
}

void SolveODE::HardcodeTaylorCoefficient(const int &coefficientNum) {
    if (coefficientNum > 0) {
        ublas::vector<double> TaylorCoefficient = ublas::vector<double>(coefficientNum, 0);
        for (int i = 0; i < coefficientNum; ++i) {
            switch (i % 4) {
                case 0:
                    TaylorCoefficient(i) = _yCauchy;
                    break;
                case 1:
                    TaylorCoefficient(i) = -_yCauchy + cos(_xCauchy);
                    break;
                case 2:
                    TaylorCoefficient(i) = _yCauchy - cos(_xCauchy) - sin(_xCauchy);
                    break;
                case 3:
                    TaylorCoefficient(i) = -_yCauchy + sin(_xCauchy);
                    break;
            }
        }
        _TaylorCoefficient = TaylorCoefficient;
    }
}

void SolveODE::BuildTaylorSolution() {
    int size = _TaylorCoefficient.size();
    Formula TaylorSolution = Formula(to_string(_yCauchy));

    for (int i = 1; i < size; ++i) {
        string powerString;
        double coefficientValue = _TaylorCoefficient(i) / factorial(i);

        if (_xCauchy >= 0) powerString = "(x - " + to_string(_xCauchy) + ")^" + to_string(i);
        else powerString = "(x + " + to_string(-_xCauchy) + ")^" + to_string(i);

        if (coefficientValue >= 0) {
            Formula coefficientFormula = to_string(coefficientValue);
            TaylorSolution = TaylorSolution + coefficientFormula * powerString;
        } else {
            Formula coefficientFormula = to_string(-coefficientValue);
            TaylorSolution = TaylorSolution - coefficientFormula * powerString;
        }
    }

    _TaylorSolution = TaylorSolution;
}

double SolveODE::EvaluateTaylorSolution(double x) {
    return _TaylorSolution.Evaluate(x);
}

void SolveODE::PrintTaylorValues() {
    int width = ceil(sqrt((double) _stepNum + 2));
    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < width; ++j) {
            int parameter = width * i + j - 2;
            if (parameter > _stepNum) break;
            cout << EvaluateTaylorSolution(_xCauchy + parameter * _step) << "\t";
        }
        cout << endl;
    }
}

void SolveODE::PrintPreciseValues(Formula &referenceFormula) {
    int width = ceil(sqrt((double) _stepNum + 2));
    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < width; ++j) {
            int parameter = width * i + j - 2;
            if (parameter > _stepNum) break;
            cout << referenceFormula.Evaluate(_xCauchy + parameter * _step) << "\t";
        }
        cout << endl;
    }
}

void SolveODE::ErrorTaylorValues(Formula &referenceFormula) {
    int width = ceil(sqrt((double) _stepNum + 2));
    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < width; ++j) {
            int parameter = width * i + j - 2;
            if (parameter > _stepNum) break;
            cout << abs(EvaluateTaylorSolution(_xCauchy + parameter * _step) -
                        referenceFormula.Evaluate(_xCauchy + parameter * _step)) << "\t";
        }
        cout << endl;
    }
}

void SolveODE::PrintValues() {
    wcout << endl << endl << L"Точка ( " << _xCauchy << " , " << _yCauchy << " )" << endl << endl;
    wcout << L"Формула: ";
    cout << _TaylorSolution.GetFormula() << endl << endl;
}

void SolveODE::InsertTaylorToAdamsMatrix() {
    for (int i = -2; i < 3; ++i) {
        double x = _xCauchy + i * _step;
        _AdamsExtrapolationMatrix(i + 2, 0) = x;
        _AdamsExtrapolationMatrix(i + 2, 1) = EvaluateTaylorSolution(x);
    }
}

void SolveODE::BuildAdamsExtrapolationMatrix() {
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5 - i; ++j) {
            if (i == 0) {
                double x = _AdamsExtrapolationMatrix(j, 0);
                double y = _AdamsExtrapolationMatrix(j, 1);
                _AdamsExtrapolationMatrix(j, 2) = _step * EvaluateDifferentialEquation(x, y);
            } else {
                _AdamsExtrapolationMatrix(j, 2 + i) =
                        _AdamsExtrapolationMatrix(j + 1, 1 + i) - _AdamsExtrapolationMatrix(j, 1 + i);
            }
        }
    }

    int stepNum = _stepNum + 2;

    for (int i = 5; i < stepNum; ++i) {
        _AdamsExtrapolationMatrix(i, 0) = _AdamsExtrapolationMatrix(i - 1, 0) + _step;

        _AdamsExtrapolationMatrix(i, 1) = _AdamsExtrapolationMatrix(i - 1, 1) + _AdamsExtrapolationMatrix(i - 1, 2) +
                                          1. / 2 * _AdamsExtrapolationMatrix(i - 2, 3) +
                                          5. / 12 * _AdamsExtrapolationMatrix(i - 3, 4) +
                                          3. / 8 * _AdamsExtrapolationMatrix(i - 4, 5) +
                                          251. / 720 * _AdamsExtrapolationMatrix(i - 5, 6);

        _AdamsExtrapolationMatrix(i, 2) =
                _step * EvaluateDifferentialEquation(_AdamsExtrapolationMatrix(i, 0), _AdamsExtrapolationMatrix(i, 1));

        for (int j = 3; j < 7; ++j) {
            _AdamsExtrapolationMatrix(i - 1, j) = _AdamsExtrapolationMatrix(i, j - 1) - _AdamsExtrapolationMatrix(i-1,j-1);
        }
    }
}



double SolveODE::EvaluateDifferentialEquation(double x, double y) {
    return _DifferentialEquation.Evaluate(x, y);
}

void SolveODE::PrintAdamsMatrix() {
    int stepNum = _stepNum - 2;

    for (int i = 0; i < stepNum; ++i) {
        for (int j = 0; j < 7; ++j) {
            cout << _AdamsExtrapolationMatrix(i, j) << "\t";
        }
        cout << endl;
    }

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 6-i; ++j) {
            cout << _AdamsExtrapolationMatrix(i + stepNum, j) << "\t";
        }
        cout << endl;
    }

}

ublas::vector<double> SolveODE::Runge_KuttaCoefficient(const double &x, const double &y) {
    ublas::vector<double> Runge_KuttaCF = ublas::vector<double>(4);
    Runge_KuttaCF(0) = _step* EvaluateDifferentialEquation(x, y);
    Runge_KuttaCF(1) = _step* EvaluateDifferentialEquation(x + _step/2, y + Runge_KuttaCF(0)/2);
    Runge_KuttaCF(2) = _step* EvaluateDifferentialEquation(x + _step/2, y + Runge_KuttaCF(1)/2);
    Runge_KuttaCF(3) = _step * EvaluateDifferentialEquation(x + _step, y + Runge_KuttaCF(2));
    return Runge_KuttaCF;
}

void SolveODE::PrintRunge_Kutta() {
    double x = _xCauchy, y = _yCauchy;

    cout << "x\ty" << endl;

    for (int i = 0; i < _stepNum; ++i) {
        cout << x << "\t" << y << endl;
        ublas::vector<double> Runge_KuttaCF = Runge_KuttaCoefficient(x, y);
        x += _step;
        y = y + 1./6*(Runge_KuttaCF(0) + 2*Runge_KuttaCF(1) + 2*Runge_KuttaCF(2) + Runge_KuttaCF(3));
    }
}

void SolveODE::EulerMethod0() {
    double x = _xCauchy, y = _yCauchy;

    cout << "x\ty" << endl;

    for (int i = 0; i < _stepNum; ++i) {
        cout << x << "\t" << y << endl;
        y = y + _step* EvaluateDifferentialEquation(x, y);
        x += _step;
    }
}

void SolveODE::EulerMethodI() {
    double x = _xCauchy, y = _yCauchy;

    cout << "x\ty" << endl;

    for (int i = 0; i < _stepNum; ++i) {
        double yHalf = y + _step/2* EvaluateDifferentialEquation(x, y);

        cout << x << "\t" << y << endl;
        y = y + _step* EvaluateDifferentialEquation(x + _step/2, yHalf);
        x += _step;
    }
}

void SolveODE::EulerMethodII() {
    double x = _xCauchy, y = _yCauchy;

    cout << "x\ty" << endl;

    for (int i = 0; i < _stepNum; ++i) {
        double yNew = y + _step* EvaluateDifferentialEquation(x, y);

        cout << x << "\t" << y << endl;
        y = y + _step/2*(EvaluateDifferentialEquation(x, y) - EvaluateDifferentialEquation(x + _step, yNew));
        x += _step;
    }
}

void SolveODE::PrintExtrapolationValues() {
    cout << "x\ty" << endl;
    for (int i = 0; i < _AdamsExtrapolationMatrix.size1(); ++i) {
        cout << _AdamsExtrapolationMatrix(i, 0) << "\t" << _AdamsExtrapolationMatrix(i, 1) << endl;
    }
}
