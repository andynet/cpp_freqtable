#include "../src/functions.cpp"

int main() {
    const char *variants_file = "../data/test/variants.txt";
    char **variants;
    uint number;

    load_variants(variants_file, &variants, &number);
}
