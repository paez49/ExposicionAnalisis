#Método de muller para aproximación de raíces, Autores: Santiago Zuñiga, Juan paez
from sympy import sympify, Symbol, lambdify, SympifyError
from cmath import sqrt

def secantMethod(f, x0, x1, error, iterations):
    for i in range(iterations):
        x2 = x1 - f(x1) * (x1 - x0) / float(f(x1) - f(x0))
        if abs(x2 - x1) < error:
            return i, x2
        x0, x1 = x1, x2
    return i, x2

def mullerMethod(f, x0, x1, x2, error, max_iter):
    for times in range(max_iter):
        # intermedios
        h0 = x1 - x0
        h1 = x2 - x1
        d0 = (f(x1)-f(x0))/h0
        d1 = (f(x2)-f(x1))/h1
        #coeficientes de la parabola
        a = (d1 - d0)/(h1 + h0)
        b = a * h1 + d1
        c = f(x2)
        #Denominador
        rad = sqrt(b**2 - 4 * a * c)
        if(abs(b+rad)>abs(b-rad)):
            deno = b+rad
        else:
            deno = b-rad
        #Hallar x3
        x3 = x2 + (-2 * c)/deno
        if abs(x3 - x2) < error:
            return times, x3
        x0, x1, x2 = x1, x2, x3

    return times, x3

def main():
    x = Symbol('x')
    try:
        f = sympify(input("Ingrese la función: "))
        f = lambdify(x, f)
        error = input("Error (default: 1e-5): ")
        error = 1e-5 if error == "" else float(error)
        iterations = input("Número max de iteraciones (dafault: 100): ")
        iterations = 100 if iterations == "" else int(iterations)
        x0 = float(input("x0: "))
        x1 = float(input("x1: "))
        x2 = float(input("x2: "))
        times, root = mullerMethod(f, x0, x1, x2, error, iterations)
        print(f"Raiz: {root}, con {times+1} iteraciones (Muller)")
        times, root = secantMethod(f, x0, x2, error, iterations)
        print(f"Raiz: {root}, con {times+1} iteraciones (Secante)")
    except SympifyError as e:
        print(f"Error parsing string: {e.expr}")
    except Exception as e:
        print("Error: ", e)

if __name__ == '__main__':
    main()
