'use strict'; // 開啟嚴格模式，避免常見的 JS 隱性錯誤

const canvas = document.getElementById('canvas'); // 取得 HTML 中的 <canvas> 元素
const ctx    = canvas.getContext('2d');            // 取得 2D 繪圖 context，所有畫圖操作都透過它

// ── 流體參數 ──────────────────────────────────────────────────────────────────
const N      = 128;    // 流體格子每軸的數量，格子越多越細緻但越吃效能
const ITER   = 4;      // Gauss-Seidel 迭代次數：越高壓力求解越準確，但越慢
const VISC   = 1e-6;   // 黏性係數：接近 0 表示像水（低黏性），調高會像蜂蜜
const DT     = 0.016;  // 每幀模擬的時間步長（對應 60fps 的一幀 ≈ 16ms）
const DAMP   = 0.998;  // 速度場每幀的衰減倍率：讓流體最終靜止，調低衰減更快
const FSCALE = 50;     // 滑鼠施力的放大倍率：調高手勢影響範圍更強
const FRAD   = 5;      // 滑鼠施力的半徑（單位：格子數）：調高影響範圍更大

// ── 粒子參數 ──────────────────────────────────────────────────────────────────
const NP       = 8000;  // 粒子總數量
const MAXSPD   = 15;    // 對應最紅色的速度上限（像素/幀），超過此值都顯示全紅
const SPRING_K = 0.008; // 回歸彈力強度：0 = 不回歸，0.02 = 很快彈回，0.008 = 緩慢飄回

// ── 格子記憶體配置 ────────────────────────────────────────────────────────────
// 格子實際大小是 (N+2)×(N+2)，多出來的 2 圈是邊界 ghost cell（用來處理邊界條件）
const SZ = (N + 2) * (N + 2);             // 格子陣列的總元素數
const IX = (i, j) => i + (N + 2) * j;    // 將 2D 座標 (i, j) 轉換成 1D 陣列索引

// u / u0：x 方向速度場（當前幀 / 上一幀，交替使用）
let u  = new Float32Array(SZ), u0 = new Float32Array(SZ);
// v / v0：y 方向速度場（當前幀 / 上一幀，交替使用）
let v  = new Float32Array(SZ), v0 = new Float32Array(SZ);
// p：壓力場（暫存，用於 project 步驟）；dv：散度場（暫存，用於 project 步驟）
let p  = new Float32Array(SZ), dv = new Float32Array(SZ);

// ── 流體求解器 ────────────────────────────────────────────────────────────────

// 設定邊界條件，讓速度在牆壁處正確反射
// b=1 表示處理 x 方向速度（u），b=2 表示 y 方向速度（v），b=0 表示純量（壓力）
function setBnd(b, x) {
  for (let i = 1; i <= N; i++) {
    // 左邊界（i=0）：x 速度取負（反射），y 速度保持同向（滑動）
    x[IX(0,   i)] = b === 1 ? -x[IX(1, i)] : x[IX(1, i)];
    // 右邊界（i=N+1）：同上邏輯，鏡射右側第一格
    x[IX(N+1, i)] = b === 1 ? -x[IX(N, i)] : x[IX(N, i)];
    // 上邊界（j=0）：y 速度取負（反射），x 速度保持同向
    x[IX(i,   0)] = b === 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
    // 下邊界（j=N+1）：同上邏輯
    x[IX(i, N+1)] = b === 2 ? -x[IX(i, N)] : x[IX(i, N)];
  }
  // 四個角落取相鄰兩格的平均，避免角落出現奇異值
  x[IX(0,   0)]   = 0.5 * (x[IX(1,   0)] + x[IX(0,   1)]);
  x[IX(N+1, 0)]   = 0.5 * (x[IX(N,   0)] + x[IX(N+1, 1)]);
  x[IX(0,   N+1)] = 0.5 * (x[IX(1, N+1)] + x[IX(0,   N)]);
  x[IX(N+1, N+1)] = 0.5 * (x[IX(N, N+1)] + x[IX(N+1, N)]);
}

// Gauss-Seidel 迭代求解線性方程組，用於 diffuse 和 project 兩個步驟
// b: 邊界類型；x: 要求解的結果陣列；x0: 已知的右側項；a, c: 方程式係數
function linSolve(b, x, x0, a, c) {
  const r = 1 / c; // 預先計算倒數，避免迴圈內重複除法
  for (let k = 0; k < ITER; k++) {       // 重複迭代 ITER 次，逐漸收斂到正確答案
    for (let j = 1; j <= N; j++) {
      for (let i = 1; i <= N; i++) {
        // 每個格子的新值 = (原始值 + a × 四鄰格的加總) / c
        // 這是在求解 (1 + 4a)·x[i,j] - a·(四鄰格) = x0[i,j] 這個方程式
        x[IX(i, j)] = (x0[IX(i, j)] + a * (
          x[IX(i-1, j)] + x[IX(i+1, j)] + // 左鄰格 + 右鄰格
          x[IX(i, j-1)] + x[IX(i, j+1)]   // 上鄰格 + 下鄰格
        )) * r;
      }
    }
    setBnd(b, x); // 每次迭代後更新邊界，確保邊界條件始終成立
  }
}

// 擴散步驟：讓速度在鄰近格子間傳播（模擬黏性）
// diff: 黏性係數；dt: 時間步長
function diffuse(b, x, x0, diff, dt) {
  const a = dt * diff * N * N; // 擴散強度：時間步長 × 黏性 × 格子總數（數值穩定化用）
  linSolve(b, x, x0, a, 1 + 4 * a); // 用隱式求解讓大步長也穩定（Stable Fluids 核心特性）
}

// 投影步驟：強制讓速度場變成「無散度」（不可壓縮），這讓流體看起來像真實液體
// ux / vy: 要修正的 x/y 速度；pArr: 壓力暫存；dvArr: 散度暫存
function project(ux, vy, pArr, dvArr) {
  const h = 1 / N; // 每個格子的寬度（在 0~1 的正規化座標中）
  for (let j = 1; j <= N; j++) {
    for (let i = 1; i <= N; i++) {
      // 計算每格的散度（速度的「擴散/收縮」程度），散度非零代表流體在壓縮或膨脹
      // 用中央差分：(右 - 左) / 2h + (下 - 上) / 2h
      dvArr[IX(i, j)] = -0.5 * h * (
        ux[IX(i+1, j)] - ux[IX(i-1, j)] + // x 方向的速度梯度
        vy[IX(i, j+1)] - vy[IX(i, j-1)]   // y 方向的速度梯度
      );
      pArr[IX(i, j)] = 0; // 壓力初始化為 0，讓 linSolve 從零開始求解
    }
  }
  setBnd(0, dvArr); setBnd(0, pArr);   // 設定散度和壓力的邊界條件（純量，b=0）
  linSolve(0, pArr, dvArr, 1, 4);      // 求解泊松方程：∇²p = 散度，得到壓力場
  for (let j = 1; j <= N; j++) {
    for (let i = 1; i <= N; i++) {
      // 用壓力梯度修正速度，讓速度場變成無散度（去掉壓縮/膨脹的分量）
      ux[IX(i, j)] -= 0.5 * (pArr[IX(i+1, j)] - pArr[IX(i-1, j)]) / h; // 修正 x 速度
      vy[IX(i, j)] -= 0.5 * (pArr[IX(i, j+1)] - pArr[IX(i, j-1)]) / h; // 修正 y 速度
    }
  }
  setBnd(1, ux); setBnd(2, vy); // 更新修正後的速度場邊界
}

// 移流步驟：讓速度場「帶著自己流動」，產生漩渦和持續的流動感
// d: 結果；d0: 來源；ux/vy: 用來移流的速度場
function advect(b, d, d0, ux, vy, dt) {
  const dt0 = dt * N; // 時間步長轉換成格子單位（N 是座標縮放因子）
  for (let j = 1; j <= N; j++) {
    for (let i = 1; i <= N; i++) {
      // 反向追蹤：現在位置 (i,j) 的流體，一幀前在哪裡？
      // 往反方向走 dt0 × 速度，就是流體的來源位置
      const x = Math.max(0.5, Math.min(N + 0.5, i - dt0 * ux[IX(i, j)])); // 來源 x（限制在有效範圍內）
      const y = Math.max(0.5, Math.min(N + 0.5, j - dt0 * vy[IX(i, j)])); // 來源 y（限制在有效範圍內）
      const i0 = Math.floor(x), i1 = i0 + 1; // 來源位置左側格子索引、右側格子索引
      const j0 = Math.floor(y), j1 = j0 + 1; // 來源位置上方格子索引、下方格子索引
      const s1 = x - i0, s0 = 1 - s1;        // x 方向雙線性插值權重（s0=左側佔比，s1=右側佔比）
      const t1 = y - j0, t0 = 1 - t1;        // y 方向雙線性插值權重（t0=上方佔比，t1=下方佔比）
      // 用雙線性插值從四個鄰格取得來源位置的值，填入當前格子
      d[IX(i, j)] =
        s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) + // 左上 × s0·t0 + 左下 × s0·t1
        s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);  // 右上 × s1·t0 + 右下 × s1·t1
    }
  }
  setBnd(b, d); // 更新結果的邊界條件
}

// 每幀執行一次完整的流體模擬步驟
function fluidStep() {
  // ── 擴散階段 ──
  // 交換 u ↔ u0（u0 變成輸入，u 變成輸出的空間），然後對 x 速度做擴散
  [u, u0] = [u0, u]; diffuse(1, u, u0, VISC, DT);
  // 同上，對 y 速度做擴散
  [v, v0] = [v0, v]; diffuse(2, v, v0, VISC, DT);
  // 擴散後先投影一次，確保速度無散度
  project(u, v, p, dv);

  // ── 移流階段 ──
  // 再次交換，讓擴散後的速度當作「移流的輸入」
  [u, u0] = [u0, u]; [v, v0] = [v0, v];
  // 用 u0/v0（移流前的速度）來移流 u（寫入移流後的 x 速度）
  advect(1, u, u0, u0, v0, DT);
  // 用 u0/v0 來移流 v（寫入移流後的 y 速度）
  advect(2, v, v0, u0, v0, DT);
  // 移流後再次投影，確保速度場仍然無散度
  project(u, v, p, dv);

  // 對整個速度場乘以衰減係數，讓流體在沒有外力時慢慢靜止
  for (let i = 0; i < SZ; i++) { u[i] *= DAMP; v[i] *= DAMP; }
}

// ── 粒子 ─────────────────────────────────────────────────────────────────────
let W = 0, H = 0;                        // 畫布寬高（在 resize 時設定）
const px    = new Float32Array(NP);      // 每顆粒子的當前 x 位置（像素）
const py    = new Float32Array(NP);      // 每顆粒子的當前 y 位置（像素）
const restX = new Float32Array(NP);      // 每顆粒子的靜止 x 位置（彈力目標點）
const restY = new Float32Array(NP);      // 每顆粒子的靜止 y 位置（彈力目標點）

// 初始化所有粒子，隨機散布在畫布上，並記錄靜止位置
function initParticles() {
  for (let i = 0; i < NP; i++) {
    px[i] = restX[i] = Math.random() * W; // 隨機 x，同時設為靜止位置
    py[i] = restY[i] = Math.random() * H; // 隨機 y，同時設為靜止位置
  }
}

// 在流體格子中對指定螢幕座標 (sx, sy) 做雙線性插值，取得該點的速度
function sampleVel(sx, sy) {
  // 將螢幕座標轉換成格子座標（1 ~ N 的範圍），並限制在有效區間內
  const gx = Math.max(0.5, Math.min(N + 0.5, 1 + sx / W * N));
  const gy = Math.max(0.5, Math.min(N + 0.5, 1 + sy / H * N));
  const i0 = Math.floor(gx), i1 = i0 + 1; // 格子座標左側和右側的整數索引
  const j0 = Math.floor(gy), j1 = j0 + 1; // 格子座標上方和下方的整數索引
  const s1 = gx - i0, s0 = 1 - s1;        // x 雙線性插值權重
  const t1 = gy - j0, t0 = 1 - t1;        // y 雙線性插值權重
  return [
    // 插值 x 速度（u）：四個鄰格加權平均
    s0 * (t0 * u[IX(i0, j0)] + t1 * u[IX(i0, j1)]) + s1 * (t0 * u[IX(i1, j0)] + t1 * u[IX(i1, j1)]),
    // 插值 y 速度（v）：四個鄰格加權平均
    s0 * (t0 * v[IX(i0, j0)] + t1 * v[IX(i0, j1)]) + s1 * (t0 * v[IX(i1, j0)] + t1 * v[IX(i1, j1)])
  ];
}

// ── 顏色 ─────────────────────────────────────────────────────────────────────
// 根據速度大小回傳顏色：慢 = 藍色（240°），快 = 紅色（0°）
function speedColor(spd) {
  const t   = Math.min(spd / MAXSPD, 1); // 將速度正規化到 0~1（超過 MAXSPD 就截斷為 1）
  const hue = (1 - t) * 240;             // t=0（靜止）→ 240° 藍，t=1（最快）→ 0° 紅
  const lit = 35 + t * 25;              // 速度越快亮度越高（35% ~ 60%），增加視覺對比
  return `hsl(${hue | 0},100%,${lit | 0}%)`; // | 0 是快速取整數（比 Math.floor 快）
}

// ── 輸入處理 ─────────────────────────────────────────────────────────────────
// 用 Map 追蹤所有活躍的指標（滑鼠用固定 key 'mouse'，觸控用 touch.identifier）
// 每個值是 { x, y, lastX, lastY }
const pointers = new Map();

function applyPointer(id, x, y) {
  const p = pointers.get(id);
  if (p) { p.lastX = p.x; p.lastY = p.y; p.x = x; p.y = y; }
}

canvas.addEventListener('mousedown', e => {
  pointers.set('mouse', { x: e.clientX, y: e.clientY, lastX: e.clientX, lastY: e.clientY });
});
canvas.addEventListener('mousemove', e => { applyPointer('mouse', e.clientX, e.clientY); });
canvas.addEventListener('mouseup',    () => { pointers.delete('mouse'); });
canvas.addEventListener('mouseleave', () => { pointers.delete('mouse'); });

canvas.addEventListener('touchstart', e => {
  e.preventDefault();
  for (const t of e.changedTouches)
    pointers.set(t.identifier, { x: t.clientX, y: t.clientY, lastX: t.clientX, lastY: t.clientY });
}, { passive: false });
canvas.addEventListener('touchmove', e => {
  e.preventDefault();
  for (const t of e.changedTouches) applyPointer(t.identifier, t.clientX, t.clientY);
}, { passive: false });
canvas.addEventListener('touchend', e => {
  e.preventDefault();
  for (const t of e.changedTouches) pointers.delete(t.identifier);
}, { passive: false });
canvas.addEventListener('touchcancel', e => {
  for (const t of e.changedTouches) pointers.delete(t.identifier);
});

// ── 主迴圈（每幀執行一次）────────────────────────────────────────────────────
function loop() {
  // ── 1. 對所有活躍指標施加力量 ──
  for (const ptr of pointers.values()) {
    const dx = ptr.x - ptr.lastX; // 本幀指標在 x 方向的位移（像素）
    const dy = ptr.y - ptr.lastY; // 本幀指標在 y 方向的位移（像素）
    const gi = Math.round(1 + ptr.x / W * N); // 指標所在格子的 i 索引
    const gj = Math.round(1 + ptr.y / H * N); // 指標所在格子的 j 索引
    for (let dj = -FRAD; dj <= FRAD; dj++) {
      for (let di = -FRAD; di <= FRAD; di++) {
        const d = Math.sqrt(di * di + dj * dj);
        if (d > FRAD) continue;
        const ii = gi + di, jj = gj + dj;
        if (ii < 1 || ii > N || jj < 1 || jj > N) continue;
        const f = 1 - d / FRAD;
        u[IX(ii, jj)] += dx / W * FSCALE * f;
        v[IX(ii, jj)] += dy / H * FSCALE * f;
      }
    }
    // 更新「上一幀位置」，供下一幀計算 dx/dy
    ptr.lastX = ptr.x;
    ptr.lastY = ptr.y;
  }

  // ── 2. 推進流體模擬一步 ──
  fluidStep();

  // ── 3. 繪製背景（半透明深藍，產生殘影/拖尾效果）──
  ctx.fillStyle = 'rgba(2, 8, 30, 0.18)'; // 每幀用低透明度藍黑色覆蓋，讓舊粒子慢慢消失
  ctx.fillRect(0, 0, W, H);               // 覆蓋整個畫布

  // ── 4. 更新並繪製每顆粒子 ──
  for (let i = 0; i < NP; i++) {
    const [ux, vy] = sampleVel(px[i], py[i]); // 在粒子當前位置取樣流體速度

    // 計算本幀的位移：流體速度（格子單位 → 像素）+ 彈力（拉向靜止位置）
    // 流體速度：格子速度 × DT × 域寬/高（換算成像素）
    // 彈力：(靜止位置 - 當前位置) × SPRING_K（越遠拉力越大）
    const vx  = ux * DT * W + (restX[i] - px[i]) * SPRING_K;
    const vy2 = vy * DT * H + (restY[i] - py[i]) * SPRING_K;

    px[i] += vx;  // 更新粒子 x 位置
    py[i] += vy2; // 更新粒子 y 位置

    // 將粒子限制在畫布範圍內（流體的 setBnd 已處理速度反射，這裡只防止位置越界）
    if (px[i] < 0) px[i] = 0;
    else if (px[i] > W) px[i] = W;
    if (py[i] < 0) py[i] = 0;
    else if (py[i] > H) py[i] = H;

    const spd = Math.sqrt(vx * vx + vy2 * vy2); // 計算本幀速度大小（像素/幀）
    ctx.fillStyle = speedColor(spd);              // 根據速度決定顏色
    ctx.fillRect(px[i] - 1, py[i] - 1, 3, 3);   // 畫出 3×3 像素的方形粒子（比 arc 快）
  }

  requestAnimationFrame(loop); // 請求瀏覽器在下一幀再呼叫 loop（維持 60fps）
}

// ── 初始化 ───────────────────────────────────────────────────────────────────
function resize() {
  W = canvas.width  = window.innerWidth;  // 更新畫布寬度為視窗寬度
  H = canvas.height = window.innerHeight; // 更新畫布高度為視窗高度
  initParticles();                         // 重新隨機散布所有粒子（視窗大小改變時重置）
}

window.addEventListener('resize', resize); // 監聽視窗大小改變
resize(); // 第一次執行，設定畫布大小並初始化粒子
loop();   // 啟動主迴圈
