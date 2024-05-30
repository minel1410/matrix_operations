#include "Matrica.h"
#include <cstdlib>
#include <ctime>


double RandomDouble(double mini, double maxi) {
    double rang = maxi - mini;
    double random = ((double) rand()) / (double) RAND_MAX;
    return mini + (random * rang);
}

int RandomInt(int mini, int maxi) {
    return rand() % (maxi - mini + 1) + mini;
}

//  Funkcija za generisanje random matrice
Matrica RandomMatrica(int brojRedova = RandomInt(1, 50), int brojKolona = RandomInt(1, 50), bool kvadratna = false) {

    Matrica m;
    if(kvadratna)
        m = Matrica(brojRedova, brojRedova);
    else
        m = Matrica(brojRedova, brojKolona);

    for(int i = 0; i < brojRedova; i++) {
        for (int j = 0; j < brojKolona; j++)
            m[i][j] = RandomDouble(-5, 5);
    }
    return m;


}
using namespace std;
//  [1 2; 3 4]^3 * (4 * [1 2 3; 4 5 6]+[7 8; 9 1; 2 4]^T) - [1 3; 7 2]^(-1)*[7 8 1; 1 2 3] * I3
// -5*[1 2; 3 7]^(-1) + ([2 3 0; 1 2 2]*[0 5.5 32; 11 4.2 3]^T) * I2*10.1
// -3.2*[3 2; 4 7]^(-2*3) + [15.1 13.2; 0.7 2]^T
// [12 34 51; 0.37 12 5.037; 2 1 3.35]*(-2) + I3^(-3)
// I3*[2 3 5; 2 0 6; 3 0.1 2.2]^T + [1 2 5; 3 4 0; 2 5 6]^(-1) - [2 2 3; 0 2 15; 3 1 0.25]^(3 * 2 + 1 - 7 + 2 - 5*2 + 17 - 2)

int main(){

    try{
        srand(time(0));

        cout<<"Testiranje konstruktora te operatora dodjele: "<<endl;
        //Matrica T(-2,2); izuzetak
        Matrica test1;
        Matrica test2(5,5);
        Matrica test3(test2);
        Matrica test4(Matrica(2,3));
        Matrica test5 = Matrica(3,2);
        Matrica test6 = test2;

        test6 += test2;
        test2 += test6;
        test4 *= test5;
        test5 /= 2;
        cout<<test1<<endl<<test2<<endl<<test3<<endl<<test4<<endl<<test5<<endl<<test6<<endl;

        cout<<"Testiranje mnozenja matrica na 10 rendom matrica: "<<endl;
        for(int i = 0; i < 5; i++) {

            Matrica A = RandomMatrica(50,15);
            Matrica B = RandomMatrica(15,50);

            cout<<A<<endl<<endl<<A*B<<endl<<endl;

        }

        for(int i = 0; i < 5; i++)
            cout<<endl;

        cout<<"Testiranje determinante na 10 rendom matrica: "<<endl;
        for(int i = 0; i < 10; i++) {

            int brojRedova = RandomInt(1, 50);
            Matrica A = RandomMatrica(brojRedova, brojRedova);
            cout<<A;
            cout<<"Determinanta: "<<A.Determinanta()<<endl<<endl<<endl;
        }

        for(int i = 0; i < 10; i++)
            cout<<endl;

        cout<<"Testiranje funkcije za invertovanje na 10 rendom matrica: "<<endl;
        for(int i = 0; i < 5; i++) {

            int brojRedova = RandomInt(1, 50);
            Matrica A = RandomMatrica(brojRedova, brojRedova);
            cout<<A<<endl;
            Matrica Ainv = A^(-1);
            cout<<Ainv<<endl;
            cout<<A*Ainv<<endl;
        }

        cout<<"Testiranje operatora '>>' "<<endl;
        Matrica A;
        cin>>A;
        cout<<A;
    }

    catch(char const poruka[]) {
        cout<<poruka;
    }
    catch(...) {}
}
