#include <iostream>
#include <vector>
#include <queue>
#include <utility>
using namespace std;

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
    {
        for (int i=0; i<M; i++)
        {
            E[U[i]].push_back(V[i]);
            Ei[U[i]].push_back(i);
            E[V[i]].push_back(U[i]);
            Ei[V[i]].push_back(i);

            R[i] = i%D;
        }

        update_dist(D);
        for (int d=0; d<D; d++)
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
                Q.pop();

                if (F[x]!=0)
                    continue;
                F[x] = 1;

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

    City city(N, M, D, K, U, V, W);

    cerr<<"Score: "<<city.calc_score()<<endl;

    for (int i=0; i<M; i++)
        cout<<(i==0?"":" ")<<city.R[i]+1;
    cout<<endl;
}

// あああああああああ
