/* Browser smoke test for generated P4 htmlwidgets.
 * Uses only Node 24 built-ins and Chrome DevTools Protocol.
 * It proves browser-side Plotly construction and one relayout interaction;
 * it is not a substitute for manual drag/selection review.
 */
const fs = require("fs");
const http = require("http");
const path = require("path");
const { spawn } = require("child_process");
const { once } = require("events");

const root = path.resolve(__dirname, "figures", "P4");
const port = Number(process.env.SC_P4_HTTP_PORT || 8765);
const debugPort = Number(process.env.SC_P4_CDP_PORT || 9229);
const chrome = process.env.CHROME_PATH ||
  "C:\\Program Files\\Google\\Chrome\\Application\\chrome.exe";
const tempRoot = require("os").tmpdir();
const userData = fs.mkdtempSync(path.join(tempRoot, "spatialcellchat-cdp-"));

function serve() {
  const server = http.createServer((req, res) => {
    const relative = decodeURIComponent((req.url || "/").split("?")[0]).replace(/^\/+/, "");
    const file = path.resolve(root, relative || "simulated/plotly_spatialDimPlot_2d.html");
    if (!file.startsWith(root + path.sep) || !fs.existsSync(file) || fs.statSync(file).isDirectory()) {
      res.writeHead(404); res.end(); return;
    }
    const ext = path.extname(file).toLowerCase();
    const type = ext === ".html" ? "text/html" :
      ext === ".js" ? "text/javascript" : ext === ".css" ? "text/css" : "application/octet-stream";
    res.writeHead(200, { "Content-Type": type });
    fs.createReadStream(file).pipe(res);
  });
  return new Promise((resolve, reject) => {
    server.once("error", reject);
    server.listen(port, "127.0.0.1", () => resolve(server));
  });
}

async function waitFor(url, timeout = 15000) {
  const started = Date.now();
  while (Date.now() - started < timeout) {
    try {
      const response = await fetch(url);
      if (response.ok) return response.json();
    } catch (_) {}
    await new Promise(resolve => setTimeout(resolve, 100));
  }
  throw new Error(`Timed out waiting for ${url}`);
}

function cdp(wsUrl) {
  const socket = new WebSocket(wsUrl);
  let next = 1;
  const pending = new Map();
  const events = new Map();
  socket.onmessage = event => {
    const message = JSON.parse(event.data);
    if (message.id && pending.has(message.id)) {
      const waiter = pending.get(message.id); pending.delete(message.id);
      if (message.error) waiter.reject(new Error(message.error.message));
      else waiter.resolve(message.result);
    } else if (message.method) {
      const listeners = events.get(message.method) || [];
      listeners.forEach(listener => listener(message.params || {}));
    }
  };
  const opened = new Promise((resolve, reject) => {
    socket.addEventListener("open", resolve, { once: true });
    socket.addEventListener("error", reject, { once: true });
  });
  return {
    async command(method, params = {}) {
      await opened;
      return new Promise((resolve, reject) => {
        const id = next++;
        pending.set(id, { resolve, reject });
        socket.send(JSON.stringify({ id, method, params }));
      });
    },
    event(method) {
      return new Promise(resolve => {
        const listeners = events.get(method) || [];
        listeners.push(resolve); events.set(method, listeners);
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

async function waitForPlotly(client) {
  const started = Date.now();
  while (Date.now() - started < 15000) {
    const ready = await evaluate(client, `(() => {
      const gd = document.querySelector('.plotly');
      return !!(gd && window.Plotly && gd._fullData && gd._fullData.length);
    })()`);
    if (ready) return;
    await new Promise(resolve => setTimeout(resolve, 100));
  }
  throw new Error("Plotly widget did not initialize in the browser");
}

const cases = [
  ["simulated", "plotly_spatialDimPlot_2d.html", "2d"],
  ["simulated", "plotly_spatialDimPlot_3d.html", "3d"],
  ["simulated", "plotly_spatialStackPlot.html", "stack"],
  ["simulated", "plotly_spatialLRpairPlot.html", "lr"],
  ["simulated", "spatialGiPlot_3d.html", "gi3d"],
  ["visium", "plotly_spatialDimPlot_2d.html", "2d"],
  ["visium", "plotly_spatialDimPlot_3d.html", "3d"],
  ["visium", "plotly_spatialStackPlot.html", "stack"],
  ["visium", "plotly_spatialLRpairPlot.html", "lr"],
  ["visium", "spatialGiPlot_3d.html", "gi3d"],
  ["slideseq", "plotly_spatialDimPlot_2d.html", "2d"],
  ["slideseq", "plotly_spatialDimPlot_3d.html", "3d"],
  ["slideseq", "plotly_spatialStackPlot.html", "stack"],
  ["slideseq", "plotly_spatialLRpairPlot.html", "lr"],
  ["slideseq", "spatialGiPlot_3d.html", "gi3d"]
];
const expectedCells = { simulated: 20, visium: 905, slideseq: 29536 };

async function main() {
  const server = await serve();
  const browser = spawn(chrome, [
    "--headless=new", "--disable-gpu", "--no-sandbox",
    `--remote-debugging-port=${debugPort}`, `--user-data-dir=${userData}`,
    "--window-size=1280,900", "about:blank"
  ], { windowsHide: true, stdio: "ignore" });
  const results = [];
  try {
    const pages = await waitFor(`http://127.0.0.1:${debugPort}/json/list`);
    const page = pages.find(item => item.type === "page");
    if (!page) throw new Error("Chrome did not expose a page target");
    const client = cdp(page.webSocketDebuggerUrl);
    await client.command("Page.enable");
    for (const [dataset, file, kind] of cases) {
      const url = `http://127.0.0.1:${port}/${dataset}/${file}`;
      await client.command("Page.navigate", { url });
      await waitForPlotly(client);
      const info = await evaluate(client, `(() => {
        const gd = document.querySelector('.plotly');
        const data = gd._fullData || [];
        const dataTraces = data.filter(trace => trace.hoverinfo !== 'none');
        const pointCount = dataTraces.reduce((total, trace) => total + (Array.isArray(trace.x) ? trace.x.length : 0), 0);
        const zPointCount = dataTraces.reduce((total, trace) => total + (Array.isArray(trace.z) ? trace.z.length : 0), 0);
        const cellIds = [...new Set(dataTraces.flatMap(trace => Array.isArray(trace.customdata) ? trace.customdata : []))];
        return {
          plotly: typeof window.Plotly,
          traceCount: data.length,
          dataTraceCount: dataTraces.length,
          pointCount,
          zPointCount,
          uniqueCellIds: cellIds.length,
          mainSvg: !!gd.querySelector('.main-svg'),
          canvas: !!gd.querySelector('canvas'),
          sceneCount: Object.keys(gd._fullLayout || {}).filter(key => /^scene\\d*$/.test(key)).length,
          colorbarCount: dataTraces.filter(trace => trace.marker && trace.marker.showscale).length,
          layerTickCount: gd._fullLayout.scene && gd._fullLayout.scene.zaxis.ticktext
            ? gd._fullLayout.scene.zaxis.ticktext.length : 0
        };
      })()`);
      const interaction = await evaluate(client, `(() => {
        const gd = document.querySelector('.plotly');
        const update = gd._fullLayout && gd._fullLayout.scene
          ? { 'scene.camera.eye': { x: 1.8, y: 1.4, z: 1.1 } }
          : { 'xaxis.autorange': true };
        return Plotly.relayout(gd, update).then(() => ({
          relayout: true,
          sceneEye: gd._fullLayout.scene && gd._fullLayout.scene.camera && gd._fullLayout.scene.camera.eye
        }));
      })()`);
      const cells = expectedCells[dataset];
      const expectedPointCount = kind === "stack" ? cells * 3 : cells;
      const hasCellCoverage = kind === "2d" || kind === "3d" || kind === "lr" || kind === "stack"
        ? info.uniqueCellIds === cells
        : info.pointCount >= cells;
      const hasExpectedGeometry = kind === "stack"
        ? info.sceneCount > 0 && info.zPointCount === expectedPointCount &&
          info.layerTickCount === 3 && info.colorbarCount === 1 && info.dataTraceCount === 1
        : kind === "3d" || kind === "gi3d"
        ? info.sceneCount > 0 && info.zPointCount >= cells
        : info.mainSvg;
      const passed = info.plotly === "object" && info.traceCount > 0 && info.pointCount >= expectedPointCount &&
        hasCellCoverage && hasExpectedGeometry && interaction.relayout;
      let screenshotPath = "";
      if (passed && kind === "stack") {
        const screenshot = await client.command("Page.captureScreenshot", {
          format: "png", captureBeyondViewport: false
        });
        screenshotPath = `${dataset}/plotly_spatialStackPlot-browser.png`;
        fs.writeFileSync(path.join(root, screenshotPath), Buffer.from(screenshot.data, "base64"));
      }
      results.push({ dataset, file, kind, status: passed ? "PASS" : "FAIL", ...info,
        browserRelayout: interaction.relayout, screenshotPath });
      console.log(`${passed ? "[PASS]" : "[FAIL]"} ${dataset}/${file} traces=${info.traceCount} points=${info.pointCount} cells=${info.uniqueCellIds} z=${info.zPointCount}`);
    }
    client.close();
  } finally {
    if (browser.exitCode === null) {
      browser.kill();
      await once(browser, "exit");
    }
    server.close();
    const resolvedTempRoot = path.resolve(tempRoot);
    const resolvedUserData = path.resolve(userData);
    if (path.dirname(resolvedUserData) !== resolvedTempRoot ||
        !path.basename(resolvedUserData).startsWith("spatialcellchat-cdp-"))
      throw new Error("Refusing to remove an unexpected Chrome profile path");
    fs.rmSync(resolvedUserData, { recursive: true, force: true });
  }
  const out = path.resolve(__dirname, "figures", "P4", "browser-summary-2026-11-10.csv");
  const header = Object.keys(results[0]).join(",");
  const csv = [header, ...results.map(row => Object.keys(results[0]).map(key => JSON.stringify(row[key] ?? "")).join(","))].join("\n") + "\n";
  fs.writeFileSync(out, csv);
  if (results.some(row => row.status !== "PASS")) process.exitCode = 1;
}

main().catch(error => { console.error(error.stack || error); process.exitCode = 1; });
