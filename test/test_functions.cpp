#include <criterion/criterion.h>
#include "../src/functions.cpp"

Test(core, empty_test) {}

Test(core, get_num_lines_returns_correct_value) {
    uint n_lines = get_n_lines("../data/test/variants.txt");
    cr_assert(n_lines == 5);
}

Test(core, load_variant_loads) {
    // Why can't I use const in the test?
    // const char *variants_file = "../data/test/variants.txt";
    char **variants;
    uint number;

    load_variants("../data/test/variants.txt", &variants, &number);
    cr_assert(strcmp(variants[0], "B.1.1.7") == 0);
    cr_assert(strcmp(variants[1], "B.1.160") == 0);
    cr_assert(strcmp(variants[2], "B.1.177") == 0);
    cr_assert(strcmp(variants[3], "B.1.258") == 0);
}
