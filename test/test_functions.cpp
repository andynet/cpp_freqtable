#include <criterion/criterion.h>
#include "../src/functions.cpp"

Test(core, empty_test) {}

Test(core, get_num_lines_returns_correct_value) {
    uint n_lines = get_n_lines("../data/test/variants.txt");
    cr_assert(n_lines == 5);
}

Test(core, load_variant_loads) {
    const char *variants_file = "../data/test/variants.txt";
    char **variants;
    uint number;

    load_variants(variants_file, &variants, &number);
    int res = strcmp(variants[0], "B.1.1.7");
    printf("%d\n", res);
}
