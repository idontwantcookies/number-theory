#include <iostream>
#include <exception>

#include <gmpxx.h>
#include <gmp.h>


mpz_class mdc_estendido(mpz_class &x, mpz_class &y, mpz_class a, mpz_class b)
{
    /* Implementa de forma iterativa o algoritmo de Euclides estendido. O
    o mdc é retornado no ponteiro g, e os coeficientes buscados retornados nos 
    ponteiros x e y, onde x é o coeficiente do maior número e y o do menor. */
    mpz_class a_, b_, resto[3], alfa[3], beta[3], quociente;
    bool trocar;

    // garantindo que a >= b
    if (a < b) { trocar = true; a_ = b; b_ = a; } 
    else { trocar = false; a_ = a; b_ = b; }

    //atribuindo valores iniciais
    alfa[0] = 1; alfa[1] = 0;
    beta[0] = 0; beta[1] = 1;
    resto[0] = a_; resto[1] = b_;
    resto[2] = a_ % b_;

    while(true) {
        quociente = resto[0] / resto[1];
        alfa[2] = alfa[0] - quociente * alfa[1];
        beta[2] = beta[0] - quociente * beta[1];
        resto[2] = alfa[2]*a_ + beta[2]*b_;

        if (resto[2] == 0) break;

        // movendo para prox linha
        for (int i = 0; i < 2; i++) {
            resto[i] = resto[i+1];
            alfa[i] = alfa[i+1];
            beta[i] = beta[i+1];
        }
    }

    if (trocar) {x = beta[1]; y = alfa[1]; } 
    else { x = alfa[1]; y = beta[1]; }
    return resto[1];
}

bool inverso_modular(mpz_class &r, mpz_class a, mpz_class n)
{
    /* usa o algoritmo de Euclides estendido para calcular o inverso modular de
    um inteiro a mod n. */

    mpz_class g, x, y;

    g = mdc_estendido(x, y, a, n);
    // if (g != 1) throw std::invalid_argument("a e n devem ser primos entre si.");
    if (g != 1) return false;
    r = x % n;
    if (r < 0) r += n;
    return true;
}

mpz_class exp_binaria(mpz_class b, mpz_class e, mpz_class n)
{
    /* calcula b^e (mod n) de forma rápida usando a decomposição do expoente e 
    em seus algarismos na base binária.  O algoritmo tem complexidade O(log e).*/

    mpz_class A = b,
              P = 1,
              E = e;

    while (E != 0) {
        if(E % 2 == 1) {
            P = (A * P) % n;
        }
        E /= 2;
        A = (A * A) % n;
    }
    return P;
}

bool primo_simples(mpz_class n)
{
    // teste inocente de primos (MUITO lento), apenas para conferir numeros pequenos
    if(n == 2 || n == -2) return true;
    if(n < 2 && n > -2) return false;
    for(int i=2; i*i<n+1; i++) {
        if(n % i == 0) return false;
    }
    return true;
}

bool primo_fermat(mpz_class n) // algoritmo de Fermat
{
    // teste determinístico se n é primo (lento, mas 100% de acerto)
    // usado para testar números pequenos e comparar com talvez_primo
    mpz_class x, y, aux, ceiling;
    
    if(n == 2 || n == -2) return true;
    if(n % 2 == 0) return false;

    mpz_sqrt(x.get_mpz_t(), n.get_mpz_t());

    if(x * x == n) return false;
    ceiling = (n + 1) / 2;
    while (true) {
        x++;
        aux = x * x - n;
        mpz_sqrt(y.get_mpz_t(), aux.get_mpz_t());
        if(x == ceiling) return true;
        if(y * y == aux) return false;
    }
}

void pre_teste_miller(mpz_class n, mpz_class &n1, unsigned int &k, mpz_class &q)
{
    // calcula n1, k e q tais que n1 = n - 1 = 2**k * q onde q é o maior ímpar possível
    k = 0;
    n1 = n - 1;
    q = n1;
    do {
        q /= 2;
        k++;
    } while(q % 2 == 0);
}

bool teste_miller(mpz_class b, mpz_class n, mpz_class n1, unsigned int k, mpz_class q)
{
    // testa na base b se um numero n e primo (teste de Miller). n1 = n-1,
    // k e q são tais que q é ímpar e 2**k*q = n - 1.
    // retorna false se o número é DEFINITIVAMENTE composto ou true se TALVEZ seja primo.
    mpz_class i, r, mdc, x, y;

    if(n == 2 || n == -2) return true;
    if(n % 2 == 0 || (n < 2 && n > -2)) return false;
    mdc = mdc_estendido(x, y, b, n);
    if(mdc == n) return true;

    r = exp_binaria(b, q, n);
    if(r == 1 || r == n1) return true;

    for(i=1; i<k; i++) {
        r = exp_binaria(r, 2, n);
        if(r == n1) return true;
    }
    return false;
}

bool primo_miller_rabin(mpz_class n, unsigned int iter, gmp_randclass &rnd)
{
    // Faz o teste de miller-rabin iter vezes usando bases aleatorias.
    // retorna false se o número é CERTAMENTE composto e true se é provavelmente primo,
    // onde a probabilidade de falso positivo primo é da ordem de 4**-iter.
    mpz_class b, n1, q;
    unsigned int k;

    if(n == 2 || n == -2) return true;
    if(-2 < n && n < 2) return false;

    pre_teste_miller(n, n1, k, q);

    for(unsigned int i=0; i<iter; i++) {
        b = rnd.get_z_range(n);
        while(b < 2) b = rnd.get_z_range(n);
        if(!teste_miller(b, n, n1, k, q)) return false;
    }
    return true;
}

mpz_class primo_aleatorio(unsigned int b, gmp_randclass &rnd)
{
    // retorna um primo aleatorio no intervalo [2, 2^b)
    // passe log=true para imprimir em clog quantas tentativas foram feitas até achar o primo
    mpz_class x;

    if(b < 1) throw std::invalid_argument("b deve ser maior ou igual a 2.");

    do {
        x = rnd.get_z_bits(b);
        x |= 1; // garantindo que seja ímpar
    } while (!primo_miller_rabin(x, 20, rnd));
    return x;
}

void gera_chaves(mpz_class &n, mpz_class &e, mpz_class &d, gmp_randclass &rnd)
{
    // gera chaves publica (n, e) e privada (n, e)
    mpz_class p, q, totiente;
    int ok;

    p = primo_aleatorio(2048, rnd);
    q = primo_aleatorio(2048, rnd);
    n = p * q;
    totiente = (p - 1) * (q - 1);
    e = 65536;
    do {
        e++;
        ok = inverso_modular(d, e, totiente);
    } while (!ok);
}

mpz_class codifica(const char* str)
{
    // recebe uma string ascii de até 500 caracteres e a codifica em um número de base 256.
    mpz_class digito = 1, r = 0;

    for(int i=0; i<500 && str[i] != '\0'; i++) {
        r += digito * str[i];
        digito *= 256;
    }
    return r;
}

char* decodifica(mpz_class n)
{
    // recebe um número em base 256 e retorna uma string de até 500 caracteres correspondente.
    char* mensagem;
    mpz_class digito;
    int i;

    mensagem = (char*) malloc(500 * sizeof(char));
    for(i=0; i<500 && n != 0; i++) {
        digito = n % 256;
        n /= 256;
        mensagem[i] = digito.get_si();
    }
    mensagem[i] = '\0';
    return mensagem;
}

mpz_class criptografa(mpz_class M, mpz_class n, mpz_class e)
{
    // retorna C = M**e % n
    return exp_binaria(M, e, n);
}

mpz_class descriptografa(mpz_class C, mpz_class n, mpz_class d)
{
    // retorna M = C**d % n
    return exp_binaria(C, d, n);
}

mpz_class gera_primo_seguro(unsigned int b, gmp_randclass& rnd)
{
    // gera um numero primo seguro p tal que p = q * 2 + 1 onde q também é primo.
    // ambos os números contêm no máximo até b bits.
    mpz_class p, q;

    do {
        do {
            q = rnd.get_z_bits(b);
            q |= 1;
            p = q * 2 + 1;
        } while (!primo_miller_rabin(p, 2, rnd) || !primo_miller_rabin(q, 2, rnd));
    } while (!primo_miller_rabin(p, 18, rnd) || !primo_miller_rabin(q, 18, rnd));
    return p;
}
