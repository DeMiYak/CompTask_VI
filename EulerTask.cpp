//
// Created by Пользователь on 19.05.2023.
//

#include"EulerTask.h"

void EulerTask(){
    Formula referenceFormula("1/2*(e^(~(x)) + cos(x) + sin(x))");
    Formula differentialEquation("~(y) + cos(x)");

    double x0 = 0, y0 = 1;
    double step;
    int stepNum;

    wcout << L"Задание 6.4" << endl
          << L"Численное решение задачи Коши для обыкновенного\nдифференциального уравнения первого порядка" << endl
          << L"Методом Эйлера (0, I, II)";
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

    wcout << L"Метод Эйлера 0: " << endl;
    ODE.EulerMethod0();
    cout << endl << endl;

    wcout << L"Метод Эйлера I: " << endl;
    ODE.EulerMethodI();
    cout << endl << endl;

    wcout << L"Метод Эйлера II: " << endl;
    ODE.EulerMethodII();
    cout << endl << endl;
}