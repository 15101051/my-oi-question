#include<algorithm>
#include<iostream>
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<map>
#include<cctype>
#include<cstring>
using namespace std;
typedef long long ll;
ll n;
ll k;
#define MAXN 10000010
#define N 10000000
#define MOD 1000000007
bool isprime[MAXN];
ll prime[MAXN],tot = 0;
ll phi[MAXN],sumphi[MAXN];
#define I inline
#define R register
I ll power(ll a,ll b)
{
	R ll res = 1;
	while(b > 0)
	{
		if(b & 1)res = 1ll * res * a % MOD;
		a = 1ll * a * a % MOD;
		b = b >> 1;
	}
	return res;
}
I ll sum(ll q,ll n)
{
	if(q == 1)return n % MOD;
	else return 1ll * (power(q,n + 1) - 1) * power(q - 1,MOD - 2) % MOD;
}
/*struct table
{
    #define MO 19260817
    int head[MO];
    ll st[4000];
    int val[4000],nxt[4000];
    int cntnum;
    I int& operator [] (ll x)
    {
        R int modx = x % MO;
        for(R int i = head[modx];i != 0;i = nxt[i])if(st[i] == x)return val[i];
        ++cntnum;
        st[cntnum] = x;val[cntnum] = 0;nxt[cntnum] = head[modx];
        head[modx] = cntnum;
        return val[cntnum];
    }
    I bool find(ll x)
    {
		R int modx = x % MO;
        for(R int i = head[modx];i != 0;i = nxt[i])if(st[i] == x)return true;
        return false;
	}
}phi_;*/
map<ll,ll> phi_;
I ll calc(ll n)
{
	if(n == 1)return 1;
	if(n <= N)return sumphi[n];
	if(phi_.find(n) != phi_.end())return phi_[n];
	R ll sum;
	if(n % 2 == 0)sum = (n / 2) % MOD * (n % MOD + 1) % MOD;
	else sum = n % MOD * ((n + 1) / 2 % MOD) % MOD;
	for(R ll l = 2,r;l <= n;l = r + 1)
	{
		r = n / (n / l);
		sum = (sum - 1ll * (r - l + 1) % MOD * calc(n / l) % MOD + MOD) % MOD;
	}
	phi_[n] = sum;
	return sum;
}
int main()
{freopen("sign10.in","r",stdin);freopen("sign10.out","w",stdout);
	scanf("%lld%d",&n,&k);
	for(R int i = 2;i <= N;++i)isprime[i] = true;
	phi[1] = sumphi[1] = 1;
	for(R int i = 2;i <= N;++i)
	{
		if(isprime[i])prime[++tot] = i,phi[i] = i - 1;
		for(R int j = 1;j <= tot && i * prime[j] <= N;++j)
		{
			R int k = i * prime[j];
			isprime[k] = false;
			if(i % prime[j] == 0){phi[k] = phi[i] * prime[j];break;}
			else phi[k] = phi[i] * phi[prime[j]];
		}
	}
	for(R int i = 1;i <= N;++i)sumphi[i] = (sumphi[i - 1] + phi[i]) % MOD;
	R ll ans = 0;
	for(R ll l = 1,r;l <= n;l = r + 1)
	{
		r = n / (n / l);
		ans = (ans + 1ll * (sum(k,r) - sum(k,l - 1) + MOD) * calc(n / l) % MOD) % MOD;
	}
	cout << ans << endl;
	return 0;
}

