/* Browser interaction checks for the communication-flow htmlwidget. */
const fs = require("fs");
const http = require("http");
const path = require("path");
const { spawn } = require("child_process");
const { once } = require("events");

const root = path.resolve(__dirname, "figures", "P3");
const port = Number(process.env.SC_FLOW_HTTP_PORT || 8767);
const debugPort = Number(process.env.SC_FLOW_CDP_PORT || 9231);
const chrome = process.env.CHROME_PATH ||
  "C:\\Program Files\\Google\\Chrome\\Application\\chrome.exe";
const tempRoot = require("os").tmpdir();
const userData = fs.mkdtempSync(path.join(tempRoot, "spatialcellchat-flow-cdp-"));
const cases = [["simulated", 20], ["visium", 905], ["slideseq", 29536]];

function serve() {
  const server = http.createServer((req, res) => {
    const relative = decodeURIComponent((req.url || "/").split("?")[0]).replace(/^\/+/, "");
    const file = path.resolve(root, relative);
    if (!file.startsWith(root + path.sep) || !fs.existsSync(file) || fs.statSync(file).isDirectory()) {
      res.writeHead(404); res.end(); return;
    }
    const ext = path.extname(file).toLowerCase();
    const type = ext === ".html" ? "text/html" : ext === ".js" ? "text/javascript" : ext === ".css" ? "text/css" : "application/octet-stream";
    res.writeHead(200, { "Content-Type": type });
    fs.createReadStream(file).pipe(res);
  });
  return new Promise((resolve, reject) => {
    server.once("error", reject);
    server.listen(port, "127.0.0.1", () => resolve(server));
  });
}

function cdp(wsUrl) {
  const socket = new WebSocket(wsUrl);
  let next = 1;
  const pending = new Map();
  const opened = new Promise((resolve, reject) => {
    socket.addEventListener("open", resolve, { once: true });
    socket.addEventListener("error", reject, { once: true });
  });
  socket.onmessage = event => {
    const message = JSON.parse(event.data);
    if (!message.id || !pending.has(message.id)) return;
    const waiter = pending.get(message.id);
    pending.delete(message.id);
    if (message.error) waiter.reject(new Error(message.error.message));
    else waiter.resolve(message.result);
  };
  return {
    async command(method, params = {}) {
      await opened;
      return new Promise((resolve, reject) => {
        const id = next++;
        pending.set(id, { resolve, reject });
        socket.send(JSON.stringify({ id, method, params }));
      });
    },
    close() { socket.close(); }
  };
}

async function evaluate(client, expression) {
  const result = await client.command("Runtime.evaluate", {
    expression, awaitPromise: true, returnByValue: true
  });
  if (result.exceptionDetails) throw new Error(result.exceptionDetails.text || "browser evaluation failed");
  return result.result.value;
}

async function waitForExpression(client, expression, timeout = 90000) {
  const started = Date.now();
  while (Date.now() - started < timeout) {
    if (await evaluate(client, expression)) return;
    await new Promise(resolve => setTimeout(resolve, 200));
  }
  throw new Error(`Timed out waiting for expression: ${expression}`);
}

async function click(client, selector) {
  const box = await evaluate(client, `(() => {
    const node = document.querySelector(${JSON.stringify(selector)});
    if (!node) return null;
    const r = node.getBoundingClientRect();
    return {x: r.left + r.width / 2, y: r.top + r.height / 2};
  })()`);
  if (!box) throw new Error(`Missing clickable element ${selector}`);
  await client.command("Input.dispatchMouseEvent", { type: "mousePressed", x: box.x, y: box.y, button: "left", clickCount: 1 });
  await client.command("Input.dispatchMouseEvent", { type: "mouseReleased", x: box.x, y: box.y, button: "left", clickCount: 1 });
}

async function setRunning(client, running) {
  const current = await evaluate(client, `document.querySelector('.sc-flow-root')?.parentElement?.__spatialCellChatFlow?.isRunning()`);
  if (current !== running) await click(client, ".sc-flow-button");
  await waitForExpression(client, `document.querySelector('.sc-flow-root')?.parentElement?.__spatialCellChatFlow?.isRunning() === ${running}`);
}

async function setRange(client, label, value) {
  return evaluate(client, `(() => {
    const control = [...document.querySelector('.sc-flow-root').querySelectorAll('label.sc-flow-control')]
      .find(item => item.querySelector('span')?.textContent === ${JSON.stringify(label)});
    if (!control) return null;
    const input = control.querySelector('input[type=range]');
    input.value = String(${JSON.stringify(value)});
    input.dispatchEvent(new Event('input', {bubbles: true}));
    return {value: input.value, output: control.querySelector('output')?.textContent};
  })()`);
}

async function state(client) {
  return evaluate(client, `(() => {
    const root = document.querySelector('.sc-flow-root');
    const widget = root?.parentElement?.__spatialCellChatFlow;
    const canvas = root?.querySelector('canvas.sc-flow-cells');
    if (!root || !widget || !canvas || !canvas.width || !canvas.height) return null;
    const image = canvas.getContext('2d').getImageData(0, 0, canvas.width, canvas.height);
    let visible = 0;
    for (let y = 0; y < canvas.height; y += 4) for (let x = 0; x < canvas.width; x += 4) {
      if (image.data[(y * canvas.width + x) * 4 + 3] > 0) visible += 1;
    }
    return {
      running: widget.isRunning(), particles: widget.particleCount(),
      snapshot: widget.particleSnapshot(), visibleCanvasPixels: visible,
      cellCount: root.querySelector('.sc-flow-status')?.textContent.match(/^(\\d+) cells/)?.[1]
    };
  })()`);
}

function moved(a, b) {
  return a.length === b.length && a.some((point, i) => Math.hypot(point[0] - b[i][0], point[1] - b[i][1]) > 1e-9);
}

async function main() {
  const server = await new Promise((resolve, reject) => {
    const value = serve();
    value.then(resolve, reject);
  });
  const browser = spawn(chrome, ["--headless=new", "--disable-gpu", "--no-sandbox",
    `--remote-debugging-port=${debugPort}`, `--user-data-dir=${userData}`,
    "--window-size=1400,1000", "about:blank"],
  { windowsHide: true, stdio: "ignore" });
  const results = [];
  let client;
  try {
    const targetsUrl = `http://127.0.0.1:${debugPort}/json/list`;
    const started = Date.now();
    let pages;
    while (Date.now() - started < 30000) {
      try {
        const response = await fetch(targetsUrl);
        if (response.ok) { pages = await response.json(); break; }
      } catch (_) {}
      await new Promise(resolve => setTimeout(resolve, 100));
    }
    const page = (pages || []).find(item => item.type === "page");
    if (!page) throw new Error("Chrome did not expose a page target");
    client = cdp(page.webSocketDebuggerUrl);
    await client.command("Page.enable");
    for (const [dataset, expectedCells] of cases) {
      const file = `${dataset}/netVisual_CommunFlowWidget.html`;
      await client.command("Page.navigate", { url: `http://127.0.0.1:${port}/${file}` });
      await waitForExpression(client, `!!document.querySelector('.sc-flow-root')?.parentElement?.__spatialCellChatFlow`);
      const initial = await state(client);
      if (!initial || Number(initial.cellCount) !== expectedCells)
        throw new Error(`${dataset}: expected ${expectedCells} cells, got ${JSON.stringify(initial)}`);
      if (!initial.running || initial.visibleCanvasPixels < 20)
        throw new Error(`${dataset}: animation/cell canvas is not active: ${JSON.stringify(initial)}`);
      await new Promise(resolve => setTimeout(resolve, 350));
      const moving = await state(client);
      const animationMoved = moved(initial.snapshot, moving.snapshot);
      if (!animationMoved) throw new Error(`${dataset}: particles did not move while running`);

      await setRunning(client, false);
      await new Promise(resolve => setTimeout(resolve, 150));
      const stopped = await state(client);
      await new Promise(resolve => setTimeout(resolve, 250));
      const stillStopped = await state(client);
      const pauseHeld = !stopped.running && !moved(stopped.snapshot, stillStopped.snapshot);
      if (!pauseHeld) throw new Error(`${dataset}: Pause did not hold particle positions`);

      const countControl = await setRange(client, "Particles", 48);
      const particleCount = await evaluate(client, `document.querySelector('.sc-flow-root')?.parentElement?.__spatialCellChatFlow?.particleCount()`);
      const opacityControl = await setRange(client, "Point opacity", 0.65);
      const zoomControl = await setRange(client, "Zoom", 1.8);
      if (!countControl || Number(particleCount) !== 48 || opacityControl?.output !== "0.65" || zoomControl?.output !== "1.8")
        throw new Error(`${dataset}: browser controls did not update: ${JSON.stringify({countControl, particleCount, opacityControl, zoomControl})}`);

      const beforeResume = await state(client);
      await setRunning(client, true);
      await new Promise(resolve => setTimeout(resolve, 350));
      const resumed = await state(client);
      const resumeMoved = moved(beforeResume.snapshot, resumed.snapshot);
      if (!resumeMoved) throw new Error(`${dataset}: Play did not resume particle motion`);

      const screenshot = await client.command("Page.captureScreenshot", { format: "png", captureBeyondViewport: false });
      fs.writeFileSync(path.join(root, dataset, "flow-widget-browser.png"), Buffer.from(screenshot.data, "base64"));
      const result = { dataset, expectedCells, renderedCells: Number(initial.cellCount),
        canvasPixels: initial.visibleCanvasPixels, animationMoved, pauseHeld,
        resumed: resumeMoved, particleCount: Number(particleCount), controls: "PASS",
        status: "PASS", screenshot: `${dataset}/flow-widget-browser.png` };
      results.push(result);
      console.log(`[PASS] ${dataset}: cells=${result.renderedCells}, canvas=${result.canvasPixels}, particles=${result.particleCount}, animation/pause/resume/controls=PASS`);
    }
  } finally {
    if (client) client.close();
    if (browser.exitCode === null) { browser.kill(); await once(browser, "exit"); }
    server.close();
    const resolvedTempRoot = path.resolve(tempRoot);
    const resolvedUserData = path.resolve(userData);
    if (path.dirname(resolvedUserData) !== resolvedTempRoot || !path.basename(resolvedUserData).startsWith("spatialcellchat-flow-cdp-"))
      throw new Error("Refusing to remove an unexpected Chrome profile path");
    fs.rmSync(resolvedUserData, { recursive: true, force: true });
  }
  if (results.length) {
    const columns = Object.keys(results[0]);
    fs.writeFileSync(path.join(root, "commun-flow-widget-browser-summary-2026-11-10.csv"),
      [columns.join(","), ...results.map(row => columns.map(key => JSON.stringify(row[key] ?? "")).join(","))].join("\n") + "\n");
  }
  if (results.length !== cases.length || results.some(row => row.status !== "PASS")) process.exitCode = 1;
}

main().catch(error => { console.error(error.stack || error); process.exitCode = 1; });
