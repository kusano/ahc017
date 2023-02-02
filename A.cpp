#include <iostream>
#include <vector>
#include <queue>
#include <utility>
#include <chrono>
#include <cmath>
#include <cstdio>
using namespace std;
using namespace std::chrono;

const long long oo = 1'000'000'000LL;

struct City
{
    int N;
    int M;
    int D;
    int K;
    vector<vector<int>> E;
    vector<vector<int>> Ei;
    vector<long long> W;
    vector<int> R;
    vector<vector<vector<long long>>> dist;
    vector<long long> F; // *N(N-1)
    vector<int> C;

    City(int N, int M, int D, int K, const vector<int> &U, const vector<int> &V, const vector<long long> &W)
        : N(N)
        , M(M)
        , D(D)
        , K(K)
        , E(N)
        , Ei(N)
        , W(W)
        , R(M)
        , dist(D+1, vector<vector<long long>>(N, vector<long long>(N)))
        , F(D+1)
        , C(D)
    {
        for (int i=0; i<M; i++)
        {
            E[U[i]].push_back(V[i]);
            Ei[U[i]].push_back(i);
            E[V[i]].push_back(U[i]);
            Ei[V[i]].push_back(i);

            R[i] = i%D;
            C[R[i]]++;
        }

        update_dist(D);
        for (int d=0; d<D; d++)
            update_dist(d);
    }

    void set(int m, int d)
    {
        if (R[m]==d)
            return;

        C[R[m]]--;
        int old = R[m];
        R[m] = d;
        C[R[m]]++;

        update_dist(old);
        update_dist(d);
    }

    void update_dist(int k)
    {
        static priority_queue<pair<long long, int>> Q;
        static vector<char> F;

        for (int p=0; p<N; p++)
        {
            for (int i=0; i<N; i++)
                dist[k][p][i] = oo;

            F.resize(N);
            for (int i=0; i<N; i++)
                F[i] = 0;

            dist[k][p][p] = 0;
            Q.push({0, p});

            while (!Q.empty())
            {
                int x = Q.top().second;
                long long qd = -Q.top().second;
                Q.pop();

                if (F[x]!=0)
                    continue;
                F[x] = 1;

                if (qd>dist[k][p][x])
                    continue;

                for (int i=0; i<(int)E[x].size(); i++)
                {
                    int e = E[x][i];
                    int ei = Ei[x][i];
                    if (R[ei]!=k)
                    {
                        long long d = dist[k][p][x]+W[ei];
                        if (d<dist[k][p][e])
                        {
                            dist[k][p][e] = d;
                            Q.push({-d, e});
                        }
                    }
                }
            }
        }

        this->F[k] = 0;
        for (int i=0; i<N; i++)
            for (int j=0; j<N; j++)
                if (i!=j)
                    this->F[k] += dist[k][i][j]-dist[D][i][j];
    }

    long long calc_score()
    {
        long long s = 0;
        for (int k=0; k<D; k++)
            s += F[k];
        long long den = D*N*(N-1);
        s = (1000*s+den/2)/den;
        return s;
    }
};

const double TIME = 5.0;
//const double TIME = 3600.0;

const int expN = 1024;
const double expX = 16;
int expT[expN];

void my_exp_init()
{
    for (int x=0; x<expN; x++)
    {
        double x2 = (double)-x/expN*expX;
        double e = exp(x2)*0x80000000+.5;
        if (e>=0x7fffffff)
            expT[x] = 0x7fffffff;
        else
            expT[x] = int(e);
    }
}

//  exp(t)*0x80000000;
int my_exp(double x)
{
    int x2 = int(x/-expX*expN+.5);
    if (x2<0)
        return expT[0];
    if (x2>=expN)
        return expT[expN-1];
    return expT[x2];
}

int xor64() {
    static uint64_t x = 88172645463345263ULL;
    x ^= x<<13;
    x ^= x>> 7;
    x ^= x<<17;
    return int(x&0x7fffffff);
}

int main()
{
    int N, M, D, K;
    cin>>N>>M>>D>>K;
    vector<int> U(M), V(M);
    vector<long long> W(M);
    for (int i=0; i<M; i++)
    {
        cin>>U[i]>>V[i]>>W[i];
        U[i]--;
        V[i]--;
    }
    vector<int> X(N), Y(N);
    for (int i=0; i<N; i++)
        cin>>X[i]>>Y[i];

    system_clock::time_point start = system_clock::now();

    City city(N, M, D, K, U, V, W);

    my_exp_init();

    long long score = city.calc_score();
    long long best_score = score;
    vector<int> best_R = city.R;

    double temp_inv;
    int iter;
    for (iter=0; ; iter++)
    {
        if (true)
        {
            system_clock::time_point now = system_clock::now();
            double time = chrono::duration_cast<chrono::nanoseconds>(now-start).count()*1e-9/TIME;
            if (time>1.0)
                break;
            double temp = 1e6*(1.0-time);
            temp_inv = 1./temp;
        }

        int m;
        int d;
        while (true)
        {
            m = xor64()%M;
            d = xor64()%D;

            if (city.R[m]!=d && city.C[d]<K)
                break;
        }

        int old_d = city.R[m];
        city.set(m, d);

        long long score2 = city.calc_score();

        if (score2<score ||
            //exp((score-score2)*temp_inv)*0x80000000>xor64())
            my_exp((score-score2)*temp_inv)>xor64())
        {
            score = score2;

            if (score<best_score)
            {
                //cerr<<score<<endl;
                best_score = score;
                best_R = city.R;
            }
        }
        else
        {
            city.set(m, old_d);
        }
    }

    //cerr<<"Iteration: "<<iter<<endl;
    //cerr<<"Score: "<<best_score<<endl;
    //cerr<<"Time: "<<chrono::duration_cast<chrono::nanoseconds>(system_clock::now()-start).count()*1e-9<<endl;

    fprintf(stderr, " %4d %4d %2d %3d %8d %16lld\n", N, M, D, K, iter, best_score);

    for (int i=0; i<M; i++)
        cout<<(i==0?"":" ")<<best_R[i]+1;
    cout<<endl;
}

// あああああああああ
