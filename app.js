const editor = document.getElementById('editor');
const preview = document.getElementById('preview');
const themeToggle = document.getElementById('theme-toggle');

const STORAGE_KEY = 'markdown-notebook-content';
const THEME_KEY = 'markdown-notebook-theme';

const DEFAULT_CONTENT = `# 歡迎使用 Markdown 筆記本

在左側輸入文字，右側即時顯示預覽。

## 功能示範

**粗體**、*斜體*、\`行內程式碼\`

### 程式碼區塊

\`\`\`javascript
function hello() {
  console.log('Hello, World!');
}
\`\`\`

### 清單

- 項目一
- 項目二
  - 子項目

### 引用

> 這是一段引用文字

---

內容會自動儲存到瀏覽器，重新整理不會遺失。
`;

function render() {
  preview.innerHTML = marked.parse(editor.value);
}

function save() {
  localStorage.setItem(STORAGE_KEY, editor.value);
}

function loadContent() {
  const saved = localStorage.getItem(STORAGE_KEY);
  editor.value = saved !== null ? saved : DEFAULT_CONTENT;
  render();
}

function loadTheme() {
  if (localStorage.getItem(THEME_KEY) === 'dark') {
    document.body.classList.add('dark');
    themeToggle.textContent = '切換淺色模式';
  }
}

editor.addEventListener('input', () => {
  render();
  save();
});

themeToggle.addEventListener('click', () => {
  const isDark = document.body.classList.toggle('dark');
  themeToggle.textContent = isDark ? '切換淺色模式' : '切換深色模式';
  localStorage.setItem(THEME_KEY, isDark ? 'dark' : 'light');
});

loadTheme();
loadContent();
