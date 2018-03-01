/*regressions*/

#include "koeffs.h"

double calculate_middle_time(
        thrust::detail::vector_base<double, thrust::device_malloc_allocator<double>>::iterator normal_iterator,
        thrust::detail::vector_base<double, thrust::device_malloc_allocator<double>>::iterator normalIterator, unsigned long i1);

double calculate_coefficient_D(
        thrust::detail::vector_base<double, thrust::device_malloc_allocator<double>>::iterator normal_iterator,
        thrust::detail::vector_base<double, thrust::device_malloc_allocator<double>>::iterator normalIterator,
        thrust::detail::vector_base<double, thrust::device_malloc_allocator<double>>::iterator iterator1,
        thrust::detail::vector_base<double, thrust::device_malloc_allocator<double>>::iterator iterator2,
        thrust::detail::vector_base<double, thrust::device_malloc_allocator<double>>::iterator iterator3, double d,
        double d1);

void printf_eps(string file_name, vector<double> e) {
    int i;
    ofstream fout("out/" + file_name, ios_base::out);
    for (i = 0; i < e.size(); i++) {
        fout << fabs(e[i]) << "\n";
    }
    fout.close();
}

int main() {
    bool    init = 1, // если данные считаны корректно, то продолжаем работу программы, иначе делаем =0 и выходим
            flag = 1,    // корректность чтения лексемы из файла
            ex1 = 1, ex2 = 1, ex3 = 1, ex4 = 1;
    int i            //счетчик
    ;
    unsigned long             //счетчик
            N; // размерность массивов
    string Line1, Line2, Line3, Line4, Line5, Line; // строковые переменные для считывания слова из файла
    double Di,Ai, Bi, j = 0.1, S = 0, S0 = 0, a = 1, tau = 1,
            chis = 0, znam = 0, // промежуточные значения
            A = 0; //коэффициенты в уравнении регрессии

    vector<double> q_oil;    //добыча жидкости
    vector<double> q_water;    //добыча воды
    vector<double> times;    //время работы по месяцам
    vector<double> e; // массив со значениями невязок

    fstream f("data/db.txt");
    ofstream fout("out/result.txt", ios_base::out);

    // считывание входных данных из файла
    //// открытие файла

    if (f.is_open()) {
        while (flag) // в случае ошибки чтения flag = 0
        {

            if (!getline(f, Line1, '\t')) flag = 0;
            if (!getline(f, Line2, '\t')) flag = 0;
            if (!getline(f, Line3, '\t')) flag = 0;
            if (!getline(f, Line4, '\t')) flag = 0;
            if (!getline(f, Line5)) flag = 0;

            times.insert(times.end(), stof(Line3));
            q_oil.insert(q_oil.end(), stof(Line4));
            q_water.insert(q_water.end(), stof(Line5));
        }
        f.close();
    } else {
        cout << "Database not found or corrupted" << endl;
        init = 0;
    }

    if (init) // если файл открылся и данные считаны, то работаем дальше
    {

        // для расчета данных по нефти надо удалить нулевые значения.

        //проверка на некорректные значения: если значение добычи <=0 то удаляем его и соответствующее
        //ему значение времени работы скважины
        for (i = 0; i < q_oil.size(); i++) {
            if (q_oil[i] <= 0) {
                q_oil.erase(q_oil.begin() + i);
                times.erase(times.begin() + i);
                q_water.erase(q_water.begin() + i);
                if (i != 0) i--;
            }
        }
        N = q_oil.size();

        // выделение памяти под массивы на видеокарте
        // ЗАМЕЧАНИЕ при таком объявлении деструктор уже встроен, удалять в конце не надо
        thrust::device_vector<double> dev_tim(N);//время
        thrust::device_vector<double> dev_oil(N);//добыча нефти
        thrust::device_vector<double> dev_res(N);//массив для записи результата
        thrust::device_vector<double> dev_wat(N);//добыча воды
        thrust::device_vector<double> dev_buf(N);//для промежуточных результатов

        // вычисляем логарифмы величин дебита и заполняем массивы на девайсе (видеокарте) - с префиксом dev_
        for (i = 0; i < N; i++) {
            dev_oil[i] = log(q_oil[i]);
            dev_wat[i] = log(q_water[i]);
            dev_tim[i] = times[i];
        }

        // Вычисление сред значений: суммируем и делим на кол-во членов в сумме
        double middle_lnq_oil = thrust::reduce(dev_oil.begin(), dev_oil.end()) / N;
        double middle_time = calculate_middle_time(dev_tim.begin(),dev_tim.end(), N);

        // линейная регрессия имеет вид: y=A-Dt;

        double coefficient_D = calculate_coefficient_D(dev_tim.begin(),
                                                       dev_tim.end(),
                                                       dev_res.begin(),
                                                       dev_res.end(),
                                                       dev_oil.begin(), middle_time, middle_lnq_oil);
        A = middle_lnq_oil - coefficient_D * middle_time;

        // Вывод на экран и сохранение в файл
        cout << "Модель 1 для нефти: " << "f(t) = exp(- " << coefficient_D << "* t) \n" << endl;
        fout << "Модель 1 для нефти: " << "f(t) = exp(- " << coefficient_D << "* t) \n" << endl;

        //файл с невязками
        for (i = 0; i < N; i++) {
            e.insert(e.end(), A - coefficient_D * times[i] - q_oil[i]);
        }
        printf_eps("eps_O1.txt", e);

        // Модель 1 для воды.
        double middle_lnq_water = thrust::reduce(dev_wat.begin(), dev_wat.end()) / N;

        //по МНК найдем D:
        thrust::transform(dev_tim.begin(), dev_tim.end(), dev_wat.begin(), dev_res.begin(), num(middle_time, middle_lnq_water));
        chis = thrust::reduce(dev_res.begin(), dev_res.end());
        coefficient_D = chis / znam;
        A = middle_lnq_water - coefficient_D * middle_time;

        // вывод результата
        cout << "Модель 1 для жидкости: " << "f(t) = exp(- " << coefficient_D << "* t) \n" << endl;
        fout << "Модель 1 для жидкости: " << "f(t) = exp(- " << coefficient_D << "* t) \n" << endl;

        // создание файла невязок
        for (i = 0; i < N; i++) {
            e[i] = A - coefficient_D * times[i] - q_water[i];
        }

        printf_eps("eps_W1.txt", e);

        // MODEL 2
        // все как в модели 1 для воды. В уравнении регрессии замена ln(1+t) = z, но на коэффициенты это не влияет.

        cout << "Модель 2 для жидкости: " << "f(t) = (1 + t) ^ - " << coefficient_D << "\n" << endl;
        fout << "Модель 2 для жидкости: " << "f(t) = (1 + t) ^ - " << coefficient_D << "\n" << endl;

        printf_eps("eps_W2.txt", e);

        ////// MODEL 3
        // коэффиц-ты для лин регрессии уже были получены в модели 1. Наша задача - найти третий коэф-т, минимизируя невязки.

        Ai = 0.0001;
        Bi = 0.0001;
        Di = 100;    //нач знач для суммы S0

        thrust::transform(dev_oil.begin(), dev_oil.end(), dev_buf.begin(), findlog(Bi, Di));
        thrust::transform(dev_tim.begin(), dev_tim.end(), dev_buf.begin(), dev_res.begin(), findS(Ai, Bi));
        S0 = thrust::reduce(dev_res.begin(), dev_res.end()) / (N - 2);

        // Покоординатный спуск с целью поиска оптимального значения трех неизвестных коэффициентов
        while (ex1 || ex2 || ex3) {
            Di += j;

            thrust::transform(dev_oil.begin(), dev_oil.end(), dev_buf.begin(), findlog(Bi, Di));
            thrust::transform(dev_tim.begin(), dev_tim.end(), dev_buf.begin(), dev_res.begin(), findS(Ai, Bi));
            S = thrust::reduce(dev_res.begin(), dev_res.end()) / (N - 2);

            if (S0 - S >= 0) {
                S0 = S;
            } else {
                Di -= j;
                ex1 = 0;
            }

            // вторая координата
            Ai += j;
            thrust::transform(dev_oil.begin(), dev_oil.end(), dev_buf.begin(), findlog(Bi, Di));
            thrust::transform(dev_tim.begin(), dev_tim.end(), dev_buf.begin(), dev_res.begin(), findS(Ai, Bi));
            S = thrust::reduce(dev_res.begin(), dev_res.end()) / (N - 2);

            if (S - S0 <= 0) {
                S0 = S;
            } else {
                Ai -= j;
                ex2 = 0;
            }

            //третья координата
            Bi += j;
            thrust::transform(dev_oil.begin(), dev_oil.end(), dev_buf.begin(), findlog(Bi, Di));
            thrust::transform(dev_tim.begin(), dev_tim.end(), dev_buf.begin(), dev_res.begin(), findS(Ai, Bi));
            S = thrust::reduce(dev_res.begin(), dev_res.end()) / (N - 2);
            if (S - S0 <= 0) {
                S0 = S;
            } else {
                Bi -= j;
                ex3 = 0;
            }
        }
        // вывод резульата
        cout << "Модель 3 для нефти: " << "f(t) = (1 + " << Bi * Di << "* t) ^ (-" << /*1. / */Bi << ") \n" << endl;
        fout << "Модель 3 для нефти: " << "f(t) = (1 + " << Bi * Di << "* t) ^ (-" << /*1. / */Bi << ") \n" << endl;

        // подсчет и вывод невязок
        for (i = 0; i < N; i++) {
            e[i] = A - coefficient_D * times[i] - q_oil[i];
        }
        printf_eps("eps_O3.txt", e);

        //// MODEL 4
        i = 0;//ограничение числа итераций: 2000
        Ai = 0.01;
        Bi = 0.01;
        Di = 0.1;
        ex1 = 1;
        ex2 = 1;
        ex3 = 1;
        thrust::transform(dev_tim.begin(), dev_tim.end(), dev_oil.begin(), dev_res.begin(), findS4(Ai, Bi, Di, tau, a));
        S0 = thrust::reduce(dev_res.begin(), dev_res.end()) / (N - 2);

        // покоординатный спуск
        while ((ex1 || ex2 || ex3 || ex4) && i < 2000) {
            tau++;
            thrust::transform(dev_tim.begin(), dev_tim.end(), dev_oil.begin(), dev_res.begin(),
                              findS4(Ai, Bi, Di, tau, a));
            S = thrust::reduce(dev_res.begin(), dev_res.end()) / (N - 2);

            if (S - S0 <= 0) {
                S0 = S;
            } else {
                tau--;
                ex1 = 0;
            }

            a += j;
            thrust::transform(dev_tim.begin(), dev_tim.end(), dev_oil.begin(), dev_res.begin(),
                              findS4(Ai, Bi, Di, tau, a));
            S = thrust::reduce(dev_res.begin(), dev_res.end()) / (N - 2);
            if (S - S0 <= 0) {
                S0 = S;
            } else {
                a -= j;
                ex2 = 0;
            }

            //B
            Bi += j;
            thrust::transform(dev_tim.begin(), dev_tim.end(), dev_oil.begin(), dev_res.begin(),
                              findS4(Ai, Bi, Di, tau, a));
            S = thrust::reduce(dev_res.begin(), dev_res.end()) / (N - 2);
            if (S - S0 <= 0) {
                S0 = S;
            } else {
                Bi -= j;
                ex3 = 0;
            }

            //A
            Ai += j;
            thrust::transform(dev_tim.begin(), dev_tim.end(), dev_oil.begin(), dev_res.begin(),
                              findS4(Ai, Bi, Di, tau, a));
            S = thrust::reduce(dev_res.begin(), dev_res.end()) / (N - 2);
            if (S - S0 <= 0) {
                S0 = S;
            } else {
                Ai -= j;
                ex4 = 0;
            }
            i++;

        }        // end of while loop

        // вывод результата
        cout << "Модель 4 для нефти: " << "f(t) = (1 + " << Bi * Di << "* t) ^ (-" << 1. / Bi << ") \n" << endl;
        fout << "Модель 4 для нефти: " << "f(t) = (1 + " << Bi * Di << "* t) ^ (-" << 1. / Bi << ") \n"
             << endl;    // и в файл
        //невязки
        for (i = 0; i < N; i++) {
            e[i] = A - coefficient_D * times[i] - q_oil[i];
        }
        printf_eps("eps_O4.txt", e);

        e.clear();
        fout.close();
    }

    times.clear();
    q_oil.clear();
    q_water.clear();

    return 0;
}
/*
 * double coefficient_D = calculate_coefficient_D(dev_tim.begin(),
                                                       dev_tim.end(),
                                                       dev_res.begin(),
                                                       dev_res.end(),
                                                       dev_oil.begin(), middle_time, middle_lnq_oil);
 */

double calculate_coefficient_D(
        thrust::detail::vector_base<double, thrust::device_malloc_allocator<double>>::iterator dev_tim_begin,
        thrust::detail::vector_base<double, thrust::device_malloc_allocator<double>>::iterator dev_tim_end,
        thrust::detail::vector_base<double, thrust::device_malloc_allocator<double>>::iterator dev_res_begin,
        thrust::detail::vector_base<double, thrust::device_malloc_allocator<double>>::iterator dev_res_end,
        thrust::detail::vector_base<double, thrust::device_malloc_allocator<double>>::iterator dev_oil_begin,
        double middle_time,
        double middle_lnq_oil) {
    //по методу наименьших квадратов найдем D:
    thrust::transform(dev_tim_begin, dev_tim_end, dev_res_begin, den(middle_time));
    double znam = thrust::reduce(dev_res_begin, dev_res_end);
    thrust::transform(dev_tim_begin, dev_tim_end, dev_oil_begin, dev_res_begin, num(middle_time, middle_lnq_oil));
    double chis = thrust::reduce(dev_res_begin, dev_res_end);
    return chis / znam;
}

double calculate_middle_time(
        thrust::detail::vector_base<double, thrust::device_malloc_allocator<double>>::iterator begin,
        thrust::detail::vector_base<double, thrust::device_malloc_allocator<double>>::iterator end, unsigned long N) {
        return thrust::reduce(begin, end) / N;;
}