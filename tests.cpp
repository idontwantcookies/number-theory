#include <iostream>
#include <string.h>

#include <gmpxx.h>
#include <gmp.h>

#include "algoritmos.hpp"

#define N 100
#define N_MUITO_LENTO 10
#define N_DETERMINISTICO 1000
#define SEED 0
#define BITS 2048

mpz_class erros;

gmp_randclass r1(gmp_randinit_default);

void testar_euclides_estendido()
{
    std::clog << "Testando euclides estendido...\n";

    mpz_class a, b, x, y, mdc, mdc_esperado;

    for(int i=0; i<N; i++)
    {
        a = r1.get_z_bits(BITS);
        b = r1.get_z_bits(BITS);

        mpz_gcd(mdc_esperado.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
        mdc = mdc_estendido(x, y, a, b);
        if(mdc != mdc_esperado) {
            erros++;
            std::cerr << "Erro: mdc(" << a << ',' << b << ") = " << mdc_esperado;
            std::cerr << ", mas o calculado foi " << mdc << ".\n";
        }
    }
}

void testar_inverso_modular()
{
    std::clog << "Testando inverso modular...\n";

    bool tem_inverso;
    mpz_class a, n, inverso, inverso_esperado;

    for(int i=0; i<N; i++)
    {
        a = r1.get_z_bits(BITS);
        n = r1.get_z_bits(BITS);
        if(n == 0) continue;

        mpz_invert(inverso_esperado.get_mpz_t(), a.get_mpz_t(), n.get_mpz_t());

        tem_inverso = inverso_modular(inverso, a, n);

        if(!tem_inverso) continue;

        if(inverso != inverso_esperado) {
            erros++;
            std::cerr << "a = " << a << ", n = " << n << '\n';
            std::cerr << "Erro: inverso modular esperado: " << inverso_esperado << '\n';
            std::cerr << "Erro: inverso modular calculado: " << inverso << '\n';
        }
    }   
}

void testar_exp_binaria()
{
    std::clog << "Testando exponenciacao binaria...\n";

    mpz_class b, e, n, result, result_esperado;

    for(int i=0; i<N; i++)
    {
        b = r1.get_z_bits(BITS);
        e = r1.get_z_bits(BITS);
        n = r1.get_z_bits(BITS);

        mpz_powm(result_esperado.get_mpz_t(), b.get_mpz_t(), e.get_mpz_t(), n.get_mpz_t());
        result = exp_binaria(b, e, n);

        if (result != result_esperado) {
            erros++;
            std::cerr << "b = " << b << ", e = " << e << ", n = " << n << '\n';
            std::cerr << "Erro: potencia esperada: " << result_esperado << '\n';
            std::cerr << "Erro: potencia calculada: " << result << '\n';
        }
    }
}

void testar_primo_fermat()
{
    std::clog << "Testando primo deterministico (Fermat)...\n";

    int primo[2];

    for(int i=0; i<N_DETERMINISTICO; i++) {
        primo[0] = primo_simples(i);
        primo[1] = primo_fermat(i);
        if(primo[0] != primo[1]) {
            erros++;
            std::cerr << "Erro: " << i << " foi considerado " << primo[0] << " no teste simples e " 
                << primo[1] << " no teste de Fermat.\n";
        }
    }
}

void testar_teste_miller() 
{
    std::clog << "Testando teste de Miller...\n";

    mpz_class n, b, swap, n1, q;
    unsigned int k;

    for(int i=0; i<N_DETERMINISTICO; i++) {
        // o teste determinístico impede que usemos bitagem alta, por ser muito lento
        n = r1.get_z_bits(20);
        b = r1.get_z_bits(20);

        pre_teste_miller(n, n1, k, q);
        if(!teste_miller(b, n, n1, k, q) && primo_fermat(n)) {
            erros++;
            std::cerr << "Erro: " << n << " com base " << b <<  "; miller: " << teste_miller(b, n, n1, k, q) << '\n';
        }
    }
}

void testar_miller_rabin() 
{
    std::clog << "Testando teste de miller-rabin...\n";

    mpz_class n;
    bool resultado[2];

    for(int i=0; i<N; i++) {
        n = r1.get_z_bits(BITS);
        resultado[0] = primo_miller_rabin(n, 20, r1);
        resultado[1] = bool(mpz_probab_prime_p(mpz_class(n).get_mpz_t(), 20));
        if(resultado[0] != resultado[1]) {
            erros++;
            std::cerr << "Erro: " << n << " teste Miller-Rabin esperado: " << resultado[1] << "; obtido: " << resultado[0] << '\n';
        }
    }
}

void testar_primalidade_pequena() 
{
    std::clog << "Testando primalidade pequena...\n";
    // imprime numeros primos pequenos. usa apenas para depuracao
    for(int i=0; i<N_DETERMINISTICO; i++) {
        if(primo_miller_rabin(i, 20, r1)) std::clog << i << ' ';
    }
    std::clog << '\n';
}

void testar_primo_aleatorio() 
{
    std::clog << "Testando primo aleatorio...\n";

    mpz_class r;

    for(int i=0; i<N_MUITO_LENTO; i++) {
        std::clog << i+1 << " de " << N_MUITO_LENTO << "...\n";
        r = primo_aleatorio(BITS, r1);
        if(!mpz_probab_prime_p(r.get_mpz_t(), 20)) {
            erros++;
            std::cerr << "Erro: " << r << " nao e primo.\n";
        }
    }
}

void estimar_bitagem_primo_aleatorio()
{
    mpz_class p;

    for(int bits = 2; bits < 1025; bits<<=1) {
        for(int i = 0; i < 100; i++) {
            p = primo_aleatorio(bits, r1);
        }
    }
}

void testar_gera_chaves()
{
    mpz_class n, e, d; // chaves
    mpz_class x, y;

    std::clog << "Testando gerador de chaves...\n";

    for(int i=0; i<N_MUITO_LENTO/2; i++) {
        std::clog << "Gerando chave " << i+1 << " de " << N_MUITO_LENTO/2 << "...\n";
        gera_chaves(n, e, d, r1);
        for(int j=2; j<N;j++) {
            x = exp_binaria(j, e, n);
            y = exp_binaria(x, d, n);
            if(j != y) {
                erros++;
                std::cerr << "Erro: Chave invalida com j = " << j << '\n';
                break;
            }
        }
    }
}

void testar_codifica()
{
    mpz_class M;
    const char *mensagem = "Mensagem ultrassecreta!";
    std::clog << "Mensagem: " << mensagem << '\n';
    M = codifica(mensagem);
    std::clog << "Mensagem codificada: " << M << '\n';
}

void testar_decodifica()
{
    mpz_class C("3197234175043038580494332653359837269788037275689313613");
    char *mensagem;
    mensagem = decodifica(C);
    std::clog << "Mensagem decodificada: " << mensagem << '\n';
}

void testar_criptografia_completa()
{
    const char *mensagem = "Mensagem ultrassecreta!";
    char *saida;
    mpz_class n, e, d;
    mpz_class M, C, D;
    std::clog << "Gerando chaves...\n";
    gera_chaves(n, e, d, r1);
    std::clog << "Mensagem: " << mensagem << '\n';
    M = codifica(mensagem);
    std::clog << "Mensagem em ASCII: " << M << '\n';
    C = criptografa(M, n, e);
    std::clog << "Mensagem criptografada: " << C << '\n';
    D = descriptografa(C, n, d);
    std::clog << "Mensagem desencriptada em ASCII: " << D << '\n';
    saida = decodifica(D);
    std::clog << "Mensagem decodificada: " << saida << '\n';
    if (strcmp(saida, mensagem)) {
        erros++;
        std::cerr << "Erro: Mensagem decodificada difere da mensagem original.\n";
    }

}

void testar_gerar_primo_seguro()
{
    std::cout << gera_primo_seguro(8, r1) << '\n';
}

int main()
{
    r1.seed(SEED);
    std::clog << "Testes iniciados com seed = " <<  SEED << '\n';
    testar_euclides_estendido();
    testar_inverso_modular();
    testar_exp_binaria();
    testar_primalidade_pequena();
    testar_primo_fermat();
    testar_teste_miller();
    testar_miller_rabin();
    testar_primo_aleatorio();
    testar_gera_chaves();
    testar_codifica();
    testar_decodifica();
    testar_criptografia_completa();
    testar_gerar_primo_seguro();

    // estimar_bitagem_primo_aleatorio();
    std::clog << erros << " erro(s) encontrado(s).\n";
    return 0;
}
