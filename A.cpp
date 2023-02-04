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

    // R[m] を d[i] にするときのコストの変化を返す。
    vector<long long> calc_diff(int m, vector<int> d)
    {
        static vector<long long> D0U(N), D0V(N), D1U(N), D1V(N);

        dijkstra2(R[m], U[m], V[m], &D0U, &D0V);

        int old = R[m];
        vector<long long> ret;
        for (int t: d)
        {
            int old = R[m];
            R[m] = t;
            dijkstra2(t, U[m], V[m], &D1U, &D1V);

            long long diff = 0;
            for (int p=0; p<N; p++)
            {
                diff -= max(0LL, abs(D0U[p]-D0V[p])-W[m]);
                diff += max(0LL, abs(D1U[p]-D1V[p])-W[m]);
            }

            long long den = 2*D*(N-1);
            diff = (1000*diff+den/2)/den;
            ret.push_back(diff);
        }
        R[m] = old;
        return ret;
    }

    // R[m]==k の辺を使わない、pからの各頂点への距離の合計を求める。
    long long dijkstra(int k, int p, int n, vector<long long> *D_)
    {
        static vector<long long> D(N);
        static priority_queue<long long> Q;

        for (int i=0; i<N; i++)
            D[i] = oo;

        D[p] = 0;
        Q.push(-p);

        long long s = 0;
        int c = 0;

        while (!Q.empty())
        {
            long long q = -Q.top();
            Q.pop();
            int x = q&0xffff;
            long long qd = q>>16;

            if (qd>D[x])
                continue;

            s += D[x];
            c++;
            if (c>=n)
                break;

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

        if (D_!=nullptr)
            *D_ = D;

        s += (n-c)*oo;
        return s;
    }

    // R[m]==k の辺を使わない、uとvからの各頂点への距離の合計を求める。
    void dijkstra2(int k, int u, int v, vector<long long> *U, vector<long long> *V)
    {
        vector<char> up(N);
        static priority_queue<long long> Q;

        for (int i=0; i<N; i++)
        {
            up[i] = 0;
            (*U)[i] = oo;
            (*V)[i] = oo;
        }

        (*U)[u] = 0;
        up[u] = 1;
        Q.push(-(oo<<16|u));
        (*V)[v] = 0;
        up[u] = 1;
        Q.push(-(oo<<16|v));

        int c = 0;
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
                if (R[ei]!=k)
                {
                    long long du = (*U)[x]+W[ei];
                    long long dv = (*V)[x]+W[ei];
                    if (du<(*U)[e] || dv<(*V)[e])
                    {
                        if (du<(*U)[e])
                            (*U)[e] = du;
                        if (dv<(*V)[e])
                            (*V)[e] = dv;
                        up[e] = true;
                        Q.push(-((du+dv)<<16|e));
                    }
                }
            }
        }
    }

    long long calc_score_orig()
    {
        long long s = 0;
        for (int p=0; p<N; p++)
        {
            for (int k=0; k<D; k++)
                s += dijkstra(k, p, N, nullptr);
            s -= dijkstra(D, p, N, nullptr)*D;
        }

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
        vector<int> ds;

        while (true)
        {
            m = xor64()%M;

            vector<int> T;
            for (int d=0; d<D; d++)
                if (city.R[m]!=d && city.C[d]<K)
                    T.push_back(d);

            if (T.empty())
                continue;

            vector<bool> U(T.size());
            int n = min((int)T.size(), 8);
            for (int i=0; i<n; i++)
            {
                int d;
                do
                    d = xor64()%T.size();
                while (U[d]);
                ds.push_back(T[d]);
            }
            break;
        }

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
