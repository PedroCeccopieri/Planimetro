// tamanho da tela
let width = 600;
let height = 700;
// coords da nova origem
let orx = width/8;
let ory = 7/8*height;

let minDist = 1;
let r = 500; // comprimento das barras

let [a, b] = [0,0]; // coordenadas do ponto (a,b)
let [x, y] = [0,0]; // coordenadas do ponto (x,y)
let [c, d] = [0,0]; // vetor perpendicular a (x,y) - (a,b)
let u1, v1, u2, v2;

let theta1 = 0; // ângulo barra vermelha
let theta2 = 0; // ângulo barra azul
// coordenadas dos mouses translatadas
let mouseXt;
let mouseYt;
// pontos pertencentes à curva
let curva = [];
let s = 0;
// variaveis para a animação
let moveTheta1 = false;
let moveTheta2 = false;
let [anteTheta1, proxTheta1, passoTheta1] = [0,0,0];
let [anteTheta2 ,proxTheta2, passoTheta2] = [0,0,0];

let isCalc = false;
let res = 0;

let idx = 0;

let img;

function preload() {
	img = loadImage('https://media.discordapp.net/attachments/529848986367164428/1065362854317920297/png-transparent-square-frame-miscellaneous-angle-white.png?width=663&height=663');
}

function setup() {
	createCanvas(width,height);

	normButton = createButton('normalisar');
	normButton.position(0,height);
	normButton.size(100,40);
	normButton.mousePressed(normalizeCurve);

	calcButton = createButton('calcular');
	calcButton.position(110,height);
	calcButton.size(100,40);
	calcButton.mousePressed(calculaArea);
}

function draw() {
	background(0,0,0); // pinta o fundo de preto
	translate(orx,ory); // define uma nova origem

	//image(img,0,0,r/1.2,-r/1.2)

	textSize(32);
	stroke(0,0,0);
	fill(255);
	text(res, 0, -500);

	mouseXt = mouseX - orx;
	mouseYt = mouseY - ory;

	// calculando as coords
	[a,b] = get_ab(theta1);
	[x,y] = get_xy(theta1,theta2);

	if (idx >= curva.length) {
		if (isCalc) {
			console.log(res)
			noLoop();
		}
		idx = 0;
	}

	desenhaPlano(); // desenha o plano cartesiano
	desenhaCurva(); // desenha a curva no canvas
	desenhaBarras(); // desenha as barras no canvas
	
	if (curva.length > 0) goto(curva[idx].x,curva[idx].y);

	// cria um ponto quanto a barra de espaço é pressionada
	if (keyIsDown(32)) {
		if (!(curva.includes({x:mouseXt,y:mouseYt}))){
			curva.push({x:mouseXt,y:mouseYt});
		}
	}

	// animação das barras
	if ((proxTheta1 - anteTheta1 < 0 && proxTheta1 < theta1) || (proxTheta1 - anteTheta1 > 0 && proxTheta1 > theta1)) {
		theta1 += passoTheta1;
		moveTheta1 = true;
	} else {
		theta1 = proxTheta1;
		moveTheta1 = false;
	}

	if ((proxTheta2 - anteTheta2 < 0 && proxTheta2 < theta2) || (proxTheta2 - anteTheta2 > 0 && proxTheta2 > theta2)) {
		theta2 += passoTheta2;
		moveTheta2 = true;
	} else {
		theta2 = proxTheta2;
		moveTheta2 = false;
	}
	
	if (idx < curva.length) {
		if (!moveTheta1 && !moveTheta2) {
			if (isCalc) calculaPonto();
			idx++;
		}
	}
}

function get_ab(theta) {
	return [r/2 * cos(theta), r/2 * sin(theta)]; // [a, b]
}

function get_xy(theta1, theta){
	return [get_ab(theta1)[0] + r/2 * cos(theta), get_ab(theta1)[1] + r/2 * sin(theta)]; // [x, y]
}

function desenhaPlano() {
	strokeWeight(5);
	stroke(255);
	var con = 1.2
	line(0,0,r/con,0); // horizontal
	line(0,0,0,-r/con); // vertical
}

function desenhaBarras() {
	strokeWeight(5);
	// barra vermelha
	stroke(255,0,0);
	line(0,0,a,b);
	// barra azul
	stroke(0,0,255);
	line(a,b,x,y);
	// barra verde
	var mag = 25/2;
	[c, d] = [(-(y-b))/(r/2), (x-a)/(r/2)]; // vetor perpendicular
	stroke(0,255,0);
	line(x - c * mag,y - d * mag,x + c * mag,y + d * mag);

	if (isCalc) {
		stroke(255,0,255);
		line(curva[idx].x - mag*u1/sqrt(u1*u1 + v1*v1), curva[idx].y - mag*v1/sqrt(u1*u1 + v1*v1), curva[idx].x + mag*u1/sqrt(u1*u1 + v1*v1), curva[idx].y + mag*v1/sqrt(u1*u1 + v1*v1));
		stroke(255,0,255);
		line(curva[idx].x - mag*u2/sqrt(u2*u2 + v2*v2), curva[idx].y - mag*v2/sqrt(u2*u2 + v2*v2), curva[idx].x + mag*u2/sqrt(u2*u2 + v2*v2), curva[idx].y + mag*v2/sqrt(u2*u2 + v2*v2));
	}
}

function desenhaCurva() {
	// desenho da curva
	strokeWeight(5);
	stroke(0,100,255);
	point(mouseXt,mouseYt);

	beginShape();
	for (let i = 0; i < curva.length; i++) {
		noFill();
		strokeWeight(5);
		stroke(255);
		curveVertex(curva[i].x,curva[i].y);
	}
	endShape(CLOSE);

	for (let i = 0; i < curva.length; i++) {
		stroke(255,255,0);
		point(curva[i].x,curva[i].y);
	}
}

function pQuinas(x,y) {
	let dis = dist(0,0,x,y);
	let h = sqrt((r/2)**2 - (dis/2)**2);
	return [x/2+h*y/dis, y/2-h*x/dis, x/2-h*y/dis, y/2+h*x/dis];
}

function anguloVetorEixoX(x,y) {

	if (y >= 0) return acos(x/(sqrt(x**2 + y**2)));
	
	if (y < 0) return -acos(x/(sqrt(x**2 + y**2)));
}

function goto(x,y) {

	var p,q;

	[p1,p2,p3,p4] = pQuinas(x,y)

	if (dist(a,b,p1,p2) <= dist(a,b,p3,p4)) [p,q] = [p1,p2];
	else [p,q] = [p3,p4];

	anteTheta1 = theta1;
	proxTheta1 = anguloVetorEixoX(p,q);
	passoTheta1 = (proxTheta1 - theta1);

	anteTheta2 = theta2;
	proxTheta2 = anguloVetorEixoX(x-get_ab(proxTheta1)[0],y-get_ab(proxTheta1)[1]);
	passoTheta2 = (proxTheta2 - theta2);
}

function removeDuplicates(arr) {
	var unique = [];
	var u = true;
    arr.forEach(p => {
		for (let i = 0; i < unique.length; i++) {
			e = unique[i];
			if (e.x == p.x && e.y == p.y) {
				u = false;
				break;
			}
		}
		if (u) unique.push(p);
		u = true;
	});
        return unique;
}

function normalizeCurve() {

	curva = removeDuplicates(curva)

	let newCurva = []

	for (let i = 0; i < curva.length; i++) {
		var k = i+1;
		if (i == curva.length - 1) k = 0;

		var qnt = floor(dist(curva[i].x,curva[i].y,curva[k].x,curva[k].y) / minDist);
		var dx = (curva[k].x - curva[i].x) / qnt;
		var dy = (curva[k].y - curva[i].y) / qnt;

		for (let j = 0; j <= qnt; j++) {
			newCurva.push({x: curva[i].x + dx * j, y: curva[i].y + dy * j});
		}
	}

	curva = removeDuplicates(newCurva);
	s = curva.length;
	idx = 0;
}

function calculaPonto() {
	var p0,p1,p2;

	p0 = curva[(idx-1 + s) % s];
    p1 = curva[(idx + s) % s];
    p2 = curva[(idx + 1 + s) % s];

	u1 = p2.x - p1.x;
	v1 = p2.y - p1.y;

	u2 = p1.x - p0.x;
	v2 = p1.y - p0.y;
	
	var cosTheta1 = -(c*u1+d*v1)/sqrt(u1**2 + v1**2);
	var cosTheta2 = -(c*u2+d*v2)/sqrt(u2**2 + v2**2);
	res += (cosTheta1 + cosTheta2)*dist(p1.x,p1.y,p2.x,p2.y)*(r/2)/(2*sqrt(c**2 + d**2));
}

function calculaArea() {
	if (curva.length > 0) {
		isCalc = true;
		idx = 0;
	}
}

function gauss() {
	var res = 0;

	for (let i = 0; i < curva.length; i++) {
        p = curva[(i - 1 + s) % s];
        q = curva[(i + s) % s];
        res += p.x * q.y;
        res -= p.y * q.x;
	}

	return abs(res)/2
}