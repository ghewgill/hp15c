load("sprintf-0.6.js")
load("jsmat/matrix.js");
load("hp15c.js");
load("test.js");

Timers = [];

function alert(s) {
    print(s);
}

function setTimeout(fn, delay) {
    var expire = new Date().getTime() + delay;
    var i = 0;
    while (i < Timers.length && expire > Timers[i].expire) {
        i += 1;
    }
    var t = {
        expire: expire,
        fn: fn
    };
    Timers[i] = t;
    return t;
}

function clearTimeout(t) {
    for (var i = 0; i < Timers.length; i++) {
        if (Timers[i] === t) {
            Timers.splice(i, 1);
            break;
        }
    }
}

console = {
    log: function(s) {
        //print(s);
    }
};

window = {
    console: console
};

Display = {
    clear_digit: function(i) {},
    clear_digits: function() {},
    clear_shift: function() {},
    set_complex: function(on) {},
    set_comma: function(i) {},
    set_decimal: function(i) {},
    set_digit: function(i, d) {},
    set_neg: function() {},
    set_prgm: function(on) {},
    set_shift: function(mode) {},
    set_trigmode: function(mode) {},
    set_user: function(on) {}
};

init();
start_tests();
while (Timers.length > 0) {
    var now = new Date().getTime();
    if (Timers.length > 0 && now >= Timers[0].expire) {
        var fn = Timers[0].fn;
        Timers.splice(0, 1);
        fn();
    }
}
