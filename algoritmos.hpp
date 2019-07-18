bool primo_simples(mpz_class);

bool primo_fermat(mpz_class);

mpz_class mdc_estendido(mpz_class&, mpz_class&, mpz_class, mpz_class);

bool inverso_modular(mpz_class&, mpz_class, mpz_class);

mpz_class exp_binaria(mpz_class, mpz_class, mpz_class);

void pre_teste_miller(mpz_class, mpz_class&, unsigned int&, mpz_class&);

bool teste_miller(mpz_class, mpz_class, mpz_class, unsigned int, mpz_class);

bool primo_miller_rabin(mpz_class, unsigned int, gmp_randclass&);

mpz_class primo_aleatorio(unsigned int, gmp_randclass&);

void gera_chaves(mpz_class&, mpz_class&, mpz_class&, gmp_randclass&);

mpz_class codifica(const char*);

char* decodifica(mpz_class);

mpz_class criptografa(mpz_class, mpz_class, mpz_class);

mpz_class descriptografa(mpz_class, mpz_class, mpz_class);

mpz_class gera_primo_seguro(unsigned int, gmp_randclass&);
