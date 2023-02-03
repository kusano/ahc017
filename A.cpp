#include <iostream>
#include <vector>
#include <queue>
#include <utility>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <utility>
using namespace std;
using namespace std::chrono;

const long long oo = 1'000'000'000LL;

struct City
{
    int N;
    int M;
    int D;
    int K;
    vector<int> U;
    vector<int> V;
    vector<vector<int>> E;
    vector<vector<int>> Ei;
    vector<long long> W;
    vector<int> X;
    vector<int> Y;
    vector<int> R;
    vector<int> C;

    City(int N, int M, int D, int K, const vector<int> &U, const vector<int> &V, const vector<long long> &W, const vector<int> &X, const vector<int> &Y)
        : N(N)
        , M(M)
        , D(D)
        , K(K)
        , U(U)
        , V(V)
        , E(N)
        , Ei(N)
        , W(W)
        , X(X)
        , Y(Y)
        , R(M)
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
    }

    // R[m] を d にするときのコストの変化を返す。
    long long calc_diff(int m, int d)
    {
        long long diff = 0;
        static vector<long long> T(N);

        dijkstra(R[m], U[m], &T);
        for (int i=0; i<N; i++)
            diff -= T[i];
        dijkstra(R[m], V[m], &T);
        for (int i=0; i<N; i++)
            diff -= T[i];
        dijkstra(d, U[m], &T);
        for (int i=0; i<N; i++)
            diff -= T[i];
        dijkstra(d, V[m], &T);
        for (int i=0; i<N; i++)
            diff -= T[i];

        int old = R[m];
        R[m] = d;
        dijkstra(old, U[m], &T);
        for (int i=0; i<N; i++)
            diff += T[i];
        dijkstra(old, V[m], &T);
        for (int i=0; i<N; i++)
            diff += T[i];
        dijkstra(d, U[m], &T);
        for (int i=0; i<N; i++)
            diff += T[i];
        dijkstra(d, V[m], &T);
        for (int i=0; i<N; i++)
            diff += T[i];
        R[m] = old;

        long long den = 2*D*(N-1);
        diff = (1000*diff+den/2)/den;
        return diff;
    }

    // R[m]==k の辺を使わない、pからの各頂点への距離を求める。
    void dijkstra(int k, int p, vector<long long> *dist)
    {
        static priority_queue<long long> Q;

        for (int i=0; i<N; i++)
            (*dist)[i] = oo;

        (*dist)[p] = 0;
        Q.push(-p);

        while (!Q.empty())
        {
            long long q = -Q.top();
            Q.pop();
            int x = q&0xffff;
            long long qd = q>>16;

            if (qd>(*dist)[x])
                continue;

            for (int i=0; i<(int)E[x].size(); i++)
            {
                int e = E[x][i];
                int ei = Ei[x][i];
                if (R[ei]!=k)
                {
                    long long d = qd+W[ei];
                    if (d<(*dist)[e])
                    {
                        (*dist)[e] = d;
                        Q.push(-(d<<16|e));
                    }
                }
            }
        }
    }

    long long calc_score_orig()
    {
        vector<vector<vector<long long>>> dist(D+1, vector<vector<long long>>(N, vector<long long>(N)));
        for (int k=0; k<=D; k++)
            for (int p=0; p<N; p++)
                dijkstra(k, p, &dist[k][p]);

        long long s = 0;
        for (int k=0; k<D; k++)
            for (int i=0; i<N; i++)
                for (int j=0; j<N; j++)
                    s += dist[k][i][j]-dist[D][i][j];
        long long den = D*N*(N-1);
        return (1000*s+den/2)/den;
    }
};

const double TIME = 5.0;

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
    if (x>=0.0)
        return 0x7fffffff;
    if (x<-expX)
        return 0;

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

    City city(N, M, D, K, U, V, W, X, Y);

    my_exp_init();

    long long score = 0;
    long long best_score = score;
    vector<int> best_R = city.R;

    double temp_inv;
    int iter;
    for (iter=0; ; iter++)
    {
        if (iter%256==0)
        {
            system_clock::time_point now = system_clock::now();
            double time = chrono::duration_cast<chrono::nanoseconds>(now-start).count()*1e-9/TIME;
            if (time>1.0)
                break;
            double temp = 1e3*(1.0-time);
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

        long long diff = city.calc_diff(m, d);
        if (diff<0 ||
            //exp(-diff*temp_inv)*0x80000000>xor64())
            my_exp(-diff*temp_inv)>xor64())
        {
            city.C[city.R[m]]--;
            city.R[m] = d;
            city.C[city.R[m]]++;

            score += diff;

            if (score<best_score)
            {
                best_score = score;
                best_R = city.R;
            }
        }
    }

    // 多い順に置き換える。
    vector<int> CC(D);
    for (int r: best_R)
        CC[r]++;
    vector<pair<int, int>> VV;
    for (int i=0; i<D; i++)
        VV.push_back({CC[i], i});
    sort(VV.begin(), VV.end());
    for (int i=0; i<M; i++)
        best_R[i] = VV[best_R[i]].second;

#ifdef TOPCODER_LOCAL
    //cerr<<"Time: "<<chrono::duration_cast<chrono::nanoseconds>(system_clock::now()-start).count()*1e-9<<endl;

    //for (auto v: VV)
    //    cerr<<" "<<v.first;
    //cerr<<endl;

    city.R = best_R;
    fprintf(stderr, " %4d %4d %2d %3d %8d %16lld %16lld\n", N, M, D, K, iter, best_score, city.calc_score_orig());
#endif

    for (int i=0; i<M; i++)
        cout<<(i==0?"":" ")<<best_R[i]+1;
    cout<<endl;
}
