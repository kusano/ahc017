#include <iostream>
#include <vector>
#include <queue>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <algorithm>
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
    vector<int> R;
    vector<int> C;
    vector<vector<int>> dist;

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
        , R(M)
        , C(D)
        , dist(N, vector<int>(N))
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

        for (int i=0; i<N; i++)
            for (int j=0; j<N; j++)
                dist[i][j] = (X[i]-X[j])*(X[i]-X[j])+(Y[i]-Y[j])*(Y[i]-Y[j]);
    }

    void set(int m, int d)
    {
        C[R[m]]--;
        R[m] = d;
        C[R[m]]++;
    }

    // R[m] を d[i] にするときのコストの変化を返す。
    vector<long long> calc_diff(int m, vector<int> ds)
    {
        static vector<long long> DU(N), DV(N);

        dijkstra2(R[m], U[m], V[m], &DU, &DV);
        long long score_old = 0;
        for (int p=0; p<N; p++)
            score_old += max(0LL, abs(DU[p]-DV[p])-W[m]);

        int old = R[m];
        vector<long long> diffs;
        for (int d: ds)
        {
            R[m] = d;
            dijkstra2(R[m], U[m], V[m], &DU, &DV);

            long long score = 0;
            for (int p=0; p<N; p++)
                score += max(0LL, abs(DU[p]-DV[p])-W[m]);

            long long den = 2*D*(N-1);
            long long diff = (1000*(score-score_old)+den/2)/den;
            diffs.push_back(diff);
        }
        R[m] = old;
        return diffs;
    }

    // R[m]==k の辺を使わない、pからの各頂点への距離の合計を求める。
    void dijkstra(int k, int p, int n, vector<long long> *D_)
    {
        static priority_queue<long long> Q;
        vector<long long> &D = *D_;

        for (int i=0; i<N; i++)
            D[i] = oo;

        D[p] = 0;
        Q.push(-p);

        while (!Q.empty())
        {
            long long q = -Q.top();
            Q.pop();
            int x = q&0xffff;
            long long qd = q>>16;

            if (qd>D[x])
                continue;

            for (int i=0; i<(int)E[x].size(); i++)
            {
                int e = E[x][i];
                int ei = Ei[x][i];
                if (R[ei]!=k)
                {
                    long long d = qd+W[ei];
                    if (d<D[e])
                    {
                        D[e] = d;
                        Q.push(-(d<<16|e));
                    }
                }
            }
        }
        while (!Q.empty())
            Q.pop();
    }

    // R[m]==k の辺を使わない、uとvからの各頂点への距離の合計を求める。
    void dijkstra2(int k, int u, int v, vector<long long> *U_, vector<long long> *V_)
    {
        static priority_queue<long long> Q;
        static vector<char> up(N);
        static vector<long long> U(N);
        static vector<long long> V(N);

        for (int i=0; i<N; i++)
        {
            up[i] = 0;
            U[i] = oo;
            V[i] = oo;
        }

        U[u] = 0;
        Q.push(-(0<<16|u));
        up[u] = 1;
        V[v] = 0;
        Q.push(-(0<<16|v));
        up[u] = 1;

        while (!Q.empty())
        {
            int x = -Q.top()&0xffff;
            Q.pop();

            if (up[x]==0)
                continue;
            up[x] = 0;

            for (int i=0; i<(int)E[x].size(); i++)
            {
                int e = E[x][i];
                int ei = Ei[x][i];
                if (R[ei]!=k && dist[u][e]<300*300)
                {
                    long long du = U[x]+W[ei];
                    long long dv = V[x]+W[ei];
                    if (du<U[e] || dv<V[e])
                    {
                        U[e] = min(U[e], du);
                        V[e] = min(V[e], dv);
                        Q.push(-(min(du,dv)<<16|e));
                        up[e] = 1;
                    }
                }
            }
        }

        *U_ = U;
        *V_ = V;
    }

    long long calc_score()
    {
        long long s = 0;
        vector<long long> dist(N);
        for (int p=0; p<N; p++)
        {
            for (int k=0; k<D; k++)
            {
                dijkstra(k, p, N, &dist);
                for (long long d: dist)
                    s += d;
            }
            dijkstra(D, p, N, &dist);
            for (long long d: dist)
                s -= d*D;
        }

        long long den = D*N*(N-1);
        return (1000*s+den/2)/den;
    }
};

const double TIME = 5.5;

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

    double temp;
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
            temp = 1e3*max(0.0, 0.9-time);
            temp_inv = 1./temp;
        }

        int m;
        vector<int> ds;

        while (true)
        {
            m = xor64()%M;

            for (int d=0; d<D; d++)
                if (city.R[m]!=d && city.C[d]<K)
                    ds.push_back(d);

            if (!ds.empty())
                break;
        }

        for (int i=(int)ds.size()-1; i>0; i--)
            swap(ds[i], ds[xor64()%(i+1)]);
        if (ds.size()>8)
            ds.resize(8);

        vector<long long> diffs = city.calc_diff(m, ds);
        long long diff = diffs[0];
        int d = ds[0];
        for (int i=1; i<(int)ds.size(); i++)
            if (diffs[i]<diff)
            {
                diff = diffs[i];
                d = ds[i];
            }

        if (diff<0 ||
            //exp(-diff*temp_inv)*0x80000000>xor64())
            (temp>0.0 && my_exp(-diff*temp_inv)>xor64()))
        {
            city.set(m, d);

            score += diff;

            if (score<best_score)
            {
                best_score = score;
                best_R = city.R;
            }
        }
    }

#ifdef TOPCODER_LOCAL
    //cerr<<"Time: "<<chrono::duration_cast<chrono::nanoseconds>(system_clock::now()-start).count()*1e-9<<endl;

    city.R = best_R;
    fprintf(stderr, " %4d %4d %2d %3d %8d %16lld %16lld\n", N, M, D, K, iter, best_score, city.calc_score());
#endif

    for (int i=0; i<M; i++)
        cout<<(i==0?"":" ")<<best_R[i]+1;
    cout<<endl;
}
