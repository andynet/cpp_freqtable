#include "../src/functions.cpp"

int main() {
    const char *variants_file = "../data/test/variants.txt";
    char **variants = NULL;
    uint number = 0;

    load_variants(variants_file, &variants, &number);
    dealloc_variants(&variants, &number);
}
