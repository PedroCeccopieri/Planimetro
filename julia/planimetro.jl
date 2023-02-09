start = time_ns();

struct point
    x::Float64
    y::Float64
end

curva = point[];
s = 0;

c, d = 0, 0;

resGreen = 0;
resGauss = 0;

r = 500;
minDist = 0.01;

A = 0;

function main()

    # Retangulo
    push!(curva,point(50.0,-100.0));
    push!(curva,point(50.0,-50.0));
    push!(curva,point(150.0,-50.0));
    push!(curva,point(150.0,-100.0));

    global A = 100 * 50;
    println("Retangulo:\nArea esperada: ",A,"\n");

    printPontos();
    calcular();
    printResultado();

    empty!(curva);

    # Quadrado
    push!(curva,point(50.0,-150.0));
    push!(curva,point(50.0,-50.0));
    push!(curva,point(150.0,-50.0));
    push!(curva,point(150.0,-150.0));

    println("----------------------------------------");

    global A = 100^2;
    println("Quadrado:\nArea esperada: ",A,"\n");

    printPontos();
    calcular();
    printResultado();

    empty!(curva);

    # Triangulo
    push!(curva,point(50.0,-150.0));
    push!(curva,point(50.0,-50.0));
    push!(curva,point(150.0,-50.0));

    println("----------------------------------------");

    global A = 100^2 / 2;
    println("Triangulo:\nArea esperada: ",A,"\n");

    printPontos();
    calcular();
    printResultado();

    empty!(curva);

    # Trapezio
    push!(curva,point(50.0,-150.0));
    push!(curva,point(50.0,-50.0));
    push!(curva,point(250.0,-50.0));
    push!(curva,point(150.0,-150.0));

    println("----------------------------------------");

    global A = (200 + 100) * 100 / 2;
    println("Trapezio:\nArea esperada: ",A,"\n");

    printPontos();
    calcular();
    printResultado();

    empty!(curva);

    n = 6;

    # Hexagono 1
    l = 20;
    for i in collect(0:2*pi/n:2*pi - 2*pi/n)
        push!(curva,point(75 + l*cos(i),-75 - l*sin(i)));
    end

    println("----------------------------------------");

    global A = 3 * sqrt(3) * l^2 / 2;
    println("Hexagono 1:\nArea esperada: ",A,"\n");

    printPontos();
    calcular();
    printResultado();

    empty!(curva);

    # Hexagono 2
    l = 100;
    for i in collect(0:2*pi/n:2*pi - 2*pi/n)
        push!(curva,point(75 + l*cos(i),-75 - l*sin(i)));
    end

    println("----------------------------------------");

    global A = 3 * sqrt(3) * l^2 / 2;
    println("Hexagono 2:\nArea esperada: ",A,"\n");

    printPontos();
    calcular();
    printResultado();

    empty!(curva);

end

function printPontos()

    println("Pontos:\n----------------------------------------")
    for p in curva
        println("(",round(p.x, digits = 2),",",round(p.y, digits = 2),")");
    end
    println("----------------------------------------")

end

function calcular()

    global curva = removeDuplicates(curva);
    global curva = normalizeCurve(curva);
    cGreen(curva);
    cGauss(curva);
    
end

function printResultado()

    println("Green: ", resGreen, " Erro: ", (A - resGreen));
    println("Gauss: ", resGauss, " Erro: ", (A - resGauss));

end

function getA(theta::Float64)
    return r/2 * cos(theta);
end

function getB(theta::Float64)
    return r/2 * sin(theta);
end

function getX(theta1::Float64, theta::Float64)
    return getA(theta1) + r/2 * cos(theta);
end

function getY(theta1::Float64, theta::Float64)
    return getB(theta1) + r/2 * sin(theta);
end

function pQuinaX(x::Float64, y::Float64)
	dis = dist(0,0,x,y);
	h = sqrt((r/2)^2 - (dis/2)^2);
	return x/2 + h * y/dis;
end

function pQuinaY(x::Float64, y::Float64)
	dis = dist(0,0,x,y);
	h = sqrt((r/2)^2 - (dis/2)^2);
	return y/2 - h * x/dis;
end

function anguloVetorEixoX(x::Float64, y::Float64)

	if (y >= 0) return acos(x/(sqrt(x*x + y*y)));
    end
	if (y < 0) return -acos(x/(sqrt(x*x + y*y)));
    end

    return 0;
end

function removeDuplicates(arr::Vector{point})
    unique::Vector{point} = [];
    un::Bool = true;
    for i in arr
        for j in unique
            if (i.x == j.x && i.y == j.y)
                un = false;
                break;
            end
        end

        if (un) 
            push!(unique, i);
        end
        un = true;
    end

    return unique;
end

function normalizeCurve(curva::Vector{point})
    
    curva = removeDuplicates(curva);

    newCurva::Vector{point} = [];
    global s = length(curva);

    for i in eachindex(curva)
        
        pa = curva[i];
        pp = curva[(i % s) + 1];
        qnt::Int = floor(dist(pa.x,pa.y,pp.x,pp.y) / minDist);
        dx::Float64 = (pp.x - pa.x) / qnt;
        dy::Float64 = (pp.y - pa.y) / qnt;

        println("criando pontos entre: (",round(pa.x, digits = 2),",",round(pa.y, digits = 2),") e (",round(pp.x, digits = 2),",",round(pp.y, digits = 2),")");

        for j in collect(0:qnt)
            push!(newCurva,point(pa.x + dx * j, pa.y + dy * j));
        end
    end

    curva = removeDuplicates(newCurva);
    global s = length(curva);
    println("----------------------------------------")
    println("Done. A curva esta com ", s, " pontos")
    println("----------------------------------------")

    return curva;
end

function calculaPonto(idx::Int)
    
    p0 = curva[((idx - 2 + s) % s) + 1];
    p1 = curva[idx];
    p2 = curva[((idx + s) % s) + 1];
    
    u1 = p2.x - p1.x;
    v1 = p2.y - p1.y;
    
    u2 = p1.x - p0.x;
    v2 = p1.y - p0.y;

    
    cosTheta1 = -(c*u1+d*v1)/sqrt(u1^2 + v1^2);
    cosTheta2 = -(c*u2+d*v2)/sqrt(u2^2 + v2^2);

    return (cosTheta1 + cosTheta2)*dist(p1.x,p1.y,p2.x,p2.y)*(r/2)/(2*sqrt(c^2 + d^2));

end

function cGreen(curva::Vector{point})
    
    res::Float64 = 0;

    for i in eachindex(curva)
        p = curva[i];
        theta1 = anguloVetorEixoX(pQuinaX(p.x,p.y),pQuinaY(p.x,p.y));
        theta2 = anguloVetorEixoX(p.x-getA(theta1),p.y-getB(theta1));
        a = getA(theta1);
        b = getB(theta1);
        x = getX(theta1,theta2);
        y = getY(theta1,theta2);
        global c = -(y-b)/(r/2);
        global d = (x-a)/(r/2);
        res += calculaPonto(i);
    end

    global resGreen = res;
end

function cGauss(curva)

    res = 0
    
    for i in eachindex(curva)
        p = curva[((i - 2 + s) % s) + 1]
        q = curva[i]
        res += p.x * q.y
        res -= p.y * q.x
    end

    global resGauss = abs(res)/2;

end

function dist(x1, y1, x2, y2)
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
end
  
main()

print("\nTempo de execucao: ", (time_ns() - start)/1000000000, " segundos");