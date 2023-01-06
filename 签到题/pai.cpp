#include<algorithm>
#include<iostream>
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<map>
#include<cctype>
#include<cstring>
using namespace std;
inline int rd()
{
	register int res = 0,f = 1;register char c = getchar();
	while(!isdigit(c)){if(c == '-')f = -1;c = getchar();}
	while(isdigit(c))res = (res << 1) + (res << 3) + c - '0',c = getchar();
	return res * f;
}
int n,k;
#define MAXN 1000010
#define MOD 1000000007
bool isprime[MAXN];
int prime[MAXN],tot = 0;
int phi[MAXN],sumphi[MAXN];
int powk[MAXN];
int main()
{
	freopen("sign7.in","r",stdin);
	freopen("sign7.out","w",stdout);
	scanf("%d%d",&n,&k);
	for(int i = 2;i <= n;++i)isprime[i] = true;
	phi[1] = 1;
	for(int i = 2;i <= n;++i)
	{
		if(isprime[i])
		{
			prime[++tot] = i;
			phi[i] = i - 1;
		}
		for(int j = 1;j <= tot && i * prime[j] <= n;++j)
		{
			int k = i * prime[j];
			isprime[k] = false;
			if(i % prime[j] == 0)
			{
				phi[k] = phi[i] * prime[j];
				break;
			}
			else
			{
				phi[k] = phi[i] * phi[prime[j]];
			}
		}
	}
	for(int i = 1;i <= n;++i)sumphi[i] = (sumphi[i - 1] + phi[i]) % MOD;
	powk[0] = 1;
	int ans = 0;
	for(int i = 1;i <= n;++i)
	{
		powk[i] = 1ll * powk[i - 1] * k % MOD;
		ans = (ans + 1ll * powk[i] * sumphi[n / i] % MOD) % MOD;
	}
	cout << ans << endl;
	return 0;
}

