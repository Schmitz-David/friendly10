from copy import deepcopy
import factor
l=[[[[5,2]],[[5,1]],[[3,2]]]]

def tupleToList(tuples):
    newlist=[]
    for t in tuples:
        newlist.append([t[0],t[1]])
    return newlist

def min_abundancy(prime_set):                   # prueft, ob die Abundanz der Zahl bereits groesser ist als 9/5
    p,q=(prime_set[0][-1][0]**(prime_set[0][-1][1]+1)-1)//(prime_set[0][-1][0]-1),prime_set[0][-1][0]**prime_set[0][-1][1]
    for i in prime_set[1]:
        p*=i[0]**i[1]
    for i in prime_set[2]:
        q*=i[0]**i[1]
    if(p==q):
        print(prime_set)                        # Bekannter von 10 gefunden!
    return p<=q

def iterate():                                  # fuehrt einen Schritt des Algorithmus durch.
    for i in range(len(l)-1,-1,-1):

        # entfernt zu grosse Elemente
        prod=(l[i][0][-1][0]**(l[i][0][-1][1]+1)-1)//(l[i][0][-1][0]-1)
        for a in l[i][2]:
            b=0
            while(prod%a[0]==0 and b<a[1]):
                b+=1
                prod//=a[0]
        for a in l[i][0]+l[i][1][1:]:
            prod*=a[0]**(a[1]+a[1]%2)
        if(l[i][0][0][1]>14 or prod>10**112):    
            del l[i]
            continue

        # Vielfachheit der letzten Primzahl erhoehen (Fall 1)
        new_prime_set=deepcopy(l[i])
        new_prime_set[0][-1][1]+=2
        l.append(new_prime_set)
        
        # Vielfachheit der letzten Primzahl festsetzen (Fall 2)
        p=l[i]
        p[2].append(deepcopy(p[0][-1]))

        p[2][-1][1]-=p[1][0][1]
        del p[1][0]
        # Die drei Listen werden entsprechend angepasst.
        factorization=fact(p[0][-1][0],p[0][-1][1]+1)
        q_index,r_index=0,0
        not_working=False
        for f in factorization: # Die neuen Primfaktoren werden eingefuegt.
            while(q_index<len(p[1])-1 and p[1][q_index][0]<f[0]):
                q_index+=1
            while(r_index<len(p[2])-1 and p[2][r_index][0]<f[0]):
                r_index+=1
            if(len(p[2])>0 and p[2][r_index][0]==f[0]):
                resulting_power=p[2][r_index][1]-f[1]
                if(resulting_power<0):
                    not_working=True        # Kandidat wird aus l geloescht.
                    break
                p[2][r_index][1]=resulting_power
                continue

            if(len(p[1])>0 and p[1][q_index][0]==f[0]):
                p[1][q_index][1]+=f[1]
            elif (len(p[1])>0):
                p[1].insert(q_index,f)
            else:
                p[1]=[f]
        if(not_working or not min_abundancy(l[i])):
            del l[i]
            continue
        p[0].append(deepcopy(p[1][0]))
        p[0][-1][1]+=p[0][-1][1]%2              # auf naechste gerade Zahl vergroessern

def factor_naiv(n):
    factorization=[]
    i=3
    while(i<=n**0.5):
        if(n%i==0):
            a=0            
            while(n%i==0):
                n=n//i
                a+=1
            factorization.append([i,a])
        i+=2
    if(n>1):
        factorization.append([n,1])
    return factorization

def show(n):
    #for prime_set in l:
    #    print(prime_set)
    print(len(l))
    iterate()
    if(n>1):
        show(n-1)

def poly_div(a,b):                      # teilt Polynom a durch Polynom b; funktioniert nur, wenn b|a
    c=a[-1]//b[-1]
    for i in range(len(b)):
        a[-i-1]-=c*b[-i-1]
        i=-1
    while(len(a) and a[-1]==0):
        del a[-1]
        i+=1
    if(len(a)):
        return poly_div(a,b)+i*[0]+[c]
    return [c]

cyclopolys=[[-1,1]]
def poly_cyc(n):                        # erzeugt die ersten n Kreisteilungspolynome
    for i in range(2,n+1):
        poly=[-1]+(i-1)*[0]+[1]
        for j in range(1,i):
            if(i%j==0):
                poly=poly_div(poly,cyclopolys[j-1])
        cyclopolys.append(poly)

poly_cyc(80)

def add(f,g):                           # addiert zwei Faktorisierungen
    result=[]
    i,j=1,1
    while(i<=len(f) and j<=len(g)):
        if(f[i-1][0]<g[j-1][0]):
            result+=[f[i-1]]
            i+=1
        else:
            if(f[i-1][0]>g[j-1][0]):
                result+=[g[j-1]]
                j+=1
            else:
                result+=[[f[i-1][0],f[i-1][1]+g[j-1][1]]]
                i+=1
                j+=1
    if(i<=len(f)):
        result+=f[i-1:]
    if(j<=len(g)):
        result+=g[j-1:]
    return result

def fact(p,n):
    result=[]
    for i in range(2,n+1):
        if(n%i==0):
            result=add(result,tupleToList(factor.factorize(sum(p**j*k for j,k in enumerate(cyclopolys[i-1])))))
    return result

def factor_trick17(n,mod):
    factorization=[]
    for p in factor_naiv(mod):
        if(n%p[0]==0):
            a=0            
            while(n%p[0]==0):
                n=n//p[0]
                a+=1
            factorization.append([p[0],a])
    p=2*mod+1
    while(p<=n**0.5):
        if(n%p==0):
            a=0            
            while(n%p==0):
                n=n//p
                a+=1
            factorization.append([p,a])
        p+=2*mod
    if(n>1):
        factorization.append([n,1])
    return factorization

show(40)
