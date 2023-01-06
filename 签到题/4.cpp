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
int k;
#define MAXN 3000010
#define N 3000000
bool isprime[MAXN];
int prime[MAXN],tot = 0;
int mu[MAXN],summu[MAXN];
#define MOD 1000000007
int power(int a,ll b)
{
	int res = 1;
	while(b > 0)
	{
		if(b & 1)res = 1ll * res * a % MOD;
		a = 1ll * a * a % MOD;
		b = b >> 1;
	}
	return res;
}
int sum(int q,ll n)
{
	if(q == 1)return n % MOD;
	else return 1ll * (power(q,n + 1) - 1) * power(q - 1,MOD - 2) % MOD;
}
#define I inline
#define R register
struct table
{
    #define MO 19260817
    int head[MO];
    ll st[400000];
    int val[400000],nxt[400000];
    int cntnum;
    I int& operator [] (ll x)
    {
        R int modx = x % MO;
        for(R int i = head[modx];i != 0;i = nxt[i])if(st[i] == x)return val[i];
        ++cntnum;st[cntnum] = x;val[cntnum] = 0;nxt[cntnum] = head[modx];head[modx] = cntnum;
        return val[cntnum];
    }
    I bool find(ll x)
    {
		R int modx = x % MO;
        for(R int i = head[modx];i != 0;i = nxt[i])if(st[i] == x)return true;
        return false;
	}
}mu_;
int calc(ll n)
{
	if(n <= N)return summu[n];
	if(mu_.find(n))return mu_[n];
	R int ans = 1;
	for(R ll l = 2,r;l <= n;l = r + 1)
	{
		r = n / (n / l);
		ans = (ans - (r - l + 1) % MOD * calc(n / l) % MOD + MOD) % MOD;
	}
	mu_[n] = ans;
	return ans;
}
I int g(ll n)
{
	R int res = 0;
	for(R ll l = 1,r;l <= n;l = r + 1)
	{
		r = n / (n / l);
		res = (res + (calc(r) - calc(l - 1) + MOD) * ((n / l) % MOD) % MOD * ((n / l) % MOD) % MOD) % MOD;
	}
	return res;
}
int main()
{freopen("Ç©µ½Ìâ/sign7.in","r",stdin);
	scanf("%lld%d",&n,&k);
	for(R int i = 2;i <= N;++i)isprime[i] = true;
	mu[1] = 1;
	for(R int i = 2;i <= N;++i)
	{
		if(isprime[i])prime[++tot] = i,mu[i] = -1;
		for(R int j = 1;j <= tot && i * prime[j] <= N;++j)
		{
			R int k = i * prime[j];
			isprime[k] = false;
			if(i % prime[j] == 0){mu[k] = 0;break;}
			else mu[k] = -mu[i];
		}
	}
	for(R int i = 1;i <= N;++i)summu[i] = (summu[i - 1] + (mu[i] + MOD) % MOD) % MOD;
	R int ans = 0;
	for(R ll l = 1,r;l <= n;l = r + 1)
	{
		r = n / (n / l);
		ans = (ans + 1ll * (sum(k,r) - sum(k,l - 1) + MOD) % MOD * g(n / l) % MOD) % MOD;
	}
	cout << 1ll * power(2,MOD - 2) * ((ans + sum(k,n) - 1) % MOD) % MOD << endl;
	return 0;
}
