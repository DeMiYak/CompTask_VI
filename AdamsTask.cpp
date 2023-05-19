//
// Created by Пользователь on 18.05.2023.
//

#include"AdamsTask.h"

void AdamsTask(){
    Formula referenceFormula("1/2*(e^(~(x)) + cos(x) + sin(x))");
    Formula differentialEquation("~(y) + cos(x)");

    double x0 = 0, y0 = 1;
    double step;
    int stepNum;

    wcout << L"Задание 6.2" << endl
          << L"Численное решение задачи Коши для обыкновенного\nдифференциального уравнения первого порядка" << endl
          << L"Экстраполяционным методом Адамса 4-го порядка";
    wcout << endl << L"Вариант 7." << endl << L"ОДУ первого порядка: " << endl;
    cout << "y'(x) = -y(x) + cos(x)";
    wcout << endl << L"С поставленной в точке (" << x0 << "," << y0 << L") задачей Коши";

    wcout << endl << L"Введите x0: ";
    cin >> x0;
    wcout << L"Введите y0: ";
    cin >> y0;
    wcout << L"Введите длину шага (h > 0): ";
    cin >> step;
    wcout << L"Введите число равноотстоящих шагов (N > 1): ";
    cin >> stepNum;

    SolveODE ODE(x0, y0, step, stepNum, differentialEquation);

    ODE.PrintValues();

    ODE.BuildAdamsExtrapolationMatrix();

    wcout << L"Вывод экстраполяционной матрицы Адамса 4-ого порядка." << endl << endl;
    ODE.PrintAdamsMatrix();
    cout << endl << endl;

    wcout << L"Вывод значений ЭМ Адамса 4-го порядка." << endl;
    ODE.PrintExtrapolationValues();
    cout << endl << endl;
}