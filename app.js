const canvas = document.getElementById('canvas');
const ctx = canvas.getContext('2d');

const GRAVITY = 0.35;
const DAMPING = 0.65;
const EMIT_COUNT = 25;

let particles = [];
let lastEmitX = -999;
let lastEmitY = -999;

function resize() {
  canvas.width = window.innerWidth;
  canvas.height = window.innerHeight;
}

class Particle {
  constructor(x, y, hue) {
    this.x = x;
    this.y = y;
    const angle = Math.random() * Math.PI * 2;
    const speed = 3 + Math.random() * 7;
    this.vx = Math.cos(angle) * speed;
    this.vy = Math.sin(angle) * speed - 2;
    this.life = 1;
    this.decay = 0.010 + Math.random() * 0.008;
    this.size = 3 + Math.random() * 4;
    this.hue = hue + (Math.random() - 0.5) * 40;
  }

  update() {
    this.vy += GRAVITY;
    this.vx *= 0.99;
    this.vy *= 0.99;
    this.x += this.vx;
    this.y += this.vy;

    if (this.x - this.size < 0) {
      this.x = this.size;
      this.vx *= -DAMPING;
    } else if (this.x + this.size > canvas.width) {
      this.x = canvas.width - this.size;
      this.vx *= -DAMPING;
    }

    if (this.y - this.size < 0) {
      this.y = this.size;
      this.vy *= -DAMPING;
    } else if (this.y + this.size > canvas.height) {
      this.y = canvas.height - this.size;
      this.vy *= -DAMPING;
    }

    this.life -= this.decay;
  }

  draw() {
    ctx.beginPath();
    ctx.arc(this.x, this.y, Math.max(0, this.size * this.life), 0, Math.PI * 2);
    ctx.fillStyle = `hsla(${this.hue}, 100%, 60%, ${this.life})`;
    ctx.fill();
  }
}

function emit(x, y) {
  const hue = Math.random() * 360;
  for (let i = 0; i < EMIT_COUNT; i++) {
    particles.push(new Particle(x, y, hue));
  }
}

function maybeEmitDrag(x, y) {
  const dx = x - lastEmitX;
  const dy = y - lastEmitY;
  if (dx * dx + dy * dy > 900) { // 30px distance threshold
    emit(x, y);
    lastEmitX = x;
    lastEmitY = y;
  }
}

function loop() {
  ctx.fillStyle = 'rgba(0, 0, 0, 0.15)';
  ctx.fillRect(0, 0, canvas.width, canvas.height);

  for (let i = particles.length - 1; i >= 0; i--) {
    particles[i].update();
    particles[i].draw();
    if (particles[i].life <= 0) {
      particles.splice(i, 1);
    }
  }

  requestAnimationFrame(loop);
}

canvas.addEventListener('mousedown', e => {
  emit(e.clientX, e.clientY);
  lastEmitX = e.clientX;
  lastEmitY = e.clientY;
});

canvas.addEventListener('mousemove', e => {
  if (e.buttons === 1) maybeEmitDrag(e.clientX, e.clientY);
});

canvas.addEventListener('touchstart', e => {
  e.preventDefault();
  for (const t of e.changedTouches) {
    emit(t.clientX, t.clientY);
    lastEmitX = t.clientX;
    lastEmitY = t.clientY;
  }
}, { passive: false });

canvas.addEventListener('touchmove', e => {
  e.preventDefault();
  for (const t of e.changedTouches) {
    maybeEmitDrag(t.clientX, t.clientY);
  }
}, { passive: false });

window.addEventListener('resize', resize);
resize();
loop();
