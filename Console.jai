
// ########  ########  #### ##    ## ########
// ##     ## ##     ##  ##  ###   ##    ##
// ##     ## ##     ##  ##  ####  ##    ##
// ########  ########   ##  ## ## ##    ##
// ##        ##   ##    ##  ##  ####    ##
// ##        ##    ##   ##  ##   ###    ##
// ##        ##     ## #### ##    ##    ##

#scope_module

PrintOptions :: struct {
    offset := 0;
    shift := 0;
    indent_every := 4;
    indent_character := "│";
};
print_options : PrintOptions;

print_indent :: (levels: int) {
    print_options.shift += 4 * levels;
}






ConsoleColor :: enum {
    FG_RESET;
    FG_BLACK;
    FG_RED;
    FG_GREEN;
    FG_YELLOW;
    FG_BLUE;
    FG_MAGENTA;
    FG_CYAN;
    FG_WHITE;
}

ansi_color :: (c: ConsoleColor = .FG_RESET) -> string {
    if c == {
        case .FG_RESET  ; return "\e[0m";
        case .FG_BLACK  ; return "\e[30m";
        case .FG_RED    ; return "\e[31m";
        case .FG_GREEN  ; return "\e[32m";
        case .FG_YELLOW ; return "\e[33m";
        case .FG_BLUE   ; return "\e[34m";
        case .FG_MAGENTA; return "\e[35m";
        case .FG_CYAN   ; return "\e[36m";
        case .FG_WHITE  ; return "\e[37m";
    }
    return ansi_color();
}



fg_color :: inline (txt: string, color: ConsoleColor = .FG_RESET) -> string {
    return sprint("%1%2%3", ansi_color(color), txt, ansi_color());
}


println_push :: (color: ConsoleColor = .FG_RESET) #expand {
    tip := cast(*Type_Info_Procedure)type_info(type_of(#this));
    text := sprint("% : % @ %", <<tip, type_of(#this), #location(#this));
    c_text := ifx color != .FG_RESET then fg_color(text, color) else text;
    println(c_text);
    print_indent(1);
    `defer println("└╴");
    `defer print_indent(-1);
}

println_push :: (text: string = "", args: ..Any, color: ConsoleColor = .FG_RESET) #expand {
    c_text := ifx color != .FG_RESET then fg_color(text, color) else text;
    println(c_text, ..args);
    print_indent(1);
    `defer println("└╴");
    `defer print_indent(-1);
}

println :: (text: string = "", args: ..Any, color: ConsoleColor = .FG_RESET) {

    c_text := ifx color != .FG_RESET then fg_color(text, color) else text;

    for 0..print_options.offset-1 {
        print(" ");
    }

    for 0..print_options.shift-1 {
        if it % print_options.indent_every == 0 {
            print(print_options.indent_character);
        } else
            print(" ");
    }

    print(c_text, ..args);
    print("\n");
}

println :: (arg: Any, color: ConsoleColor = .FG_RESET) {
    text := sprint("%", arg);
    c_text := ifx color != .FG_RESET then fg_color(text, color) else text;
    
    for 0..print_options.offset-1 {
        print(" ");
    }

    for 0..print_options.shift-1 {
        if it % print_options.indent_every == 0 {
            print(print_options.indent_character);
        } else
            print(" ");
    }
    print("%\n", c_text);
}



test_console :: () {
    println_push("Console", color = .FG_WHITE);
    println_push(.FG_BLUE);
}