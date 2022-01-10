#import "MathExtensions";
#import "Basic";

main :: () {
    verbose := false;
    args := get_command_line_arguments();

    if args.count > 1 {
        if args[1] == "-verbose" verbose = true;
    }

    run_all_tests(verbose);
}
