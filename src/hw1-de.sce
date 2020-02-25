/*
    IVP:
    y' = -y + 4x*x + 12*x + 98
    
    Initial condition:    y(0) = 1
    Interval of interest: [0; 4]

    Exact solution: 
    y = 4x*x+4x+94-93exp(-x)
*/

function ydot=F(x, y)
    ydot=-y+4*x^2+12*x+98
endfunction

function y = solve_exact(y0,x)
    C = (y0-(4*x(1)^2+4*x(1)+94))*exp(x(1))
    y = 4*x^2+4*x+94+C*exp(-x)
endfunction

function y=solve_built_in(y0,x,f)
    y = ode(y0, x(1),x, f)
endfunction

function y=solve_euler(y0,x,h,f)
    y = zeros(x) // (1xN)    
    d = zeros(size(x)(2)-1)
    y(1) = y0
    for i = 1:(size(x)(2)-1)
        d(i) = h * f(x(i), y(i))
        y(i+1) = y(i) + d(i)
    end
endfunction

function y=solve_euler_improved(y0,x,h,f)
    y = zeros(x) // (1xN)
    y(1) = y0
    for i = 1:(size(x)(2)-1)
        m1 = f(x(i), y(i))
        m2 = f(x(i + 1), y(i) + h * m1)
        y(i + 1) = (y(i) + h * (m1 + m2) / 2)
    end
endfunction

function y=solve_runge_kutta(y0,x,h,f)
    y = zeros(x)    
    y(1) = y0 
    for i = 1:(size(x)(2)-1)
        k1 = f(x(i), y(i))
        k2 = f(x(i) + h / 2, y(i) + h / 2 * k1)
        k3 = f(x(i) + h / 2, y(i) + h / 2 * k2)
        k4 = f(x(i) + h, y(i) + h * k3)

        y(i + 1) = (y(i) + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4))
    end
endfunction

function err = loc_err(exact, numerical)
    err = zeros(exact)
    for i = 1:size(err)(2)
        err(i) = abs(exact(i)-numerical(i))    
    end 
endfunction

function er = glob_err(loc_err)
    er = max(loc_err)
endfunction

/*
    Plots
*/
function plot_sol(x, values,labels,styles,wind)
    show_window(wind)
    title('Solutions of the IVP: y = -y + 4x*x + 12*x + 98', "fontsize", 2)
    legends(labels, styles, opt="lr")
    xlabel('X in Interval', "fontsize", 1)
    ylabel('Y', "fontsize", 1)
    plot2d(x, values, styles)
endfunction

function plot_loc_err_log(x, values,labels,styles,wind)
    show_window(wind)
    title('Local Error Plots: ', "fontsize", 2)
    legends(labels, styles, opt="lr")
    xlabel('X in Interval', "fontsize", 1)
    ylabel('Local Error log10', "fontsize", 1)
    plot2d(x, log10(values), styles)
endfunction


function plot_glob_err_log(x, values,labels,styles,wind)
    show_window(wind)
    title("Global Error Plots: ", "fontsize", 2)
    legends(labels, styles, opt="lr")
    xlabel('N - number of steps', "fontsize", 1)
    ylabel('Global Error log10', "fontsize", 1)
    plot2d(x, log10(values), styles)
endfunction

/*
    Solver
*/
function run_sols()
    // configurations
    I = [0,4]
    N = 50
    h = (I(2) - I(1)) / N
    x = I(1):h:I(2)
    y0=0

    // solutions
    disp("N:",N)
    exact = solve_exact(y0,x)    
    standard = solve_built_in(y0,x,F)
    standard_e = loc_err(exact, standard)
    standard_g = glob_err(standard_e')
    disp("Standard ODE global error:")
    disp(standard_g)
        
    euler = solve_euler(y0,x, h, F)
    euler_e = loc_err(exact, euler)
    euler_g = glob_err(euler_e')
    disp("Euler global error:")
    disp(euler_g)    
    
    euler_improved = solve_euler_improved(y0,x, h, F)
    euler_improved_e = loc_err(exact, euler_improved)
    euler_improved_g = glob_err(euler_improved_e')
    disp("Euler Improved global error:")
    disp(euler_improved_g)    
    
    runge_kutta = solve_runge_kutta(y0,x, h, F)
    runge_kutta_e = loc_err(exact, runge_kutta)
    runge_kutta_g = glob_err(runge_kutta_e')
    disp("Runge Kutta global error:")
    disp(runge_kutta_g)    
    
    // plottings
    sols = [exact',standard',euler',euler_improved',runge_kutta']
    loc_errs = [standard_e',euler_e',euler_improved_e',runge_kutta_e']
    labels = ['exact'; 'standard';'euler';'euler improved';'runge kutta']
    labels_e = ['standard error';'euler error';'euler improved error'; ..
    'runge kutta error']
    styles = [1,2,3,4,5]
    plot_sol(x,sols,labels,styles,1)
    plot_loc_err_log(x,loc_errs,labels_e,styles,2)
    
endfunction

function run_errs()
    // configurations
    I = [0,4]
    N = 50
    y0=0
 
    euler_errs = zeros(1,N)
    euler_improved_errs = zeros(1,N)
    runge_kutta_errs = zeros(1,N)
    for i= 1:N
        h = (I(2) - I(1)) / i
        x = I(1):h:I(2)
        
        exact = solve_exact(y0,x)    
        sol = solve_euler(y0,x,h,F)
        sol_e = loc_err(exact, sol)
        sol_g = glob_err(sol_e')
        euler_errs(i) = sol_g

        sol = solve_euler_improved(y0,x,h,F)
        sol_e = loc_err(exact, sol)
        sol_g = glob_err(sol_e')
        euler_improved_errs(i) = sol_g

        sol = solve_runge_kutta(y0,x,h,F)
        sol_e = loc_err(exact, sol)
        sol_g = glob_err(sol_e')
        runge_kutta_errs(i) = sol_g
    end
    glob_errs = [euler_errs',euler_improved_errs',runge_kutta_errs']
    labels_e = ['euler error';'euler improved error';'runge kutta error']
    styles = [1,2,3]
    plot_glob_err_log(1:N,glob_errs,labels_e,styles,3)
endfunction


run_sols()
run_errs()

