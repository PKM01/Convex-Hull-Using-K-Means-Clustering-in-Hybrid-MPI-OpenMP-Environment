#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>

using namespace std;

int main()
{
    int n;
    cout << "Enter the number of points to generate: ";
    cin >> n;

    ofstream outfile("points_" + to_string(n) + ".txt");
    if (!outfile.is_open())
    {
        cout << "Error: could not open file." << endl;
        return 1;
    }

    srand(time(NULL));
    for (int i = 0; i < n; i++)
    {
        int x = rand() % 100;
        int y = rand() % 100;
        outfile << x << " " << y << endl;
    }

    outfile.close();
    cout << "Points written to file 'points.txt'." << endl;
    return 0;
}
