/* Browser smoke for Visium and full Slide-seq Shiny selection paths. */
const fs = require("fs");
const path = require("path");
const { spawn } = require("child_process");
const { once } = require("events");

const dataset = process.env.SC_SHINY_DATASET;
const expectedCells = { visium: 905, slideseq: 29536 }[dataset];
if (!expectedCells) throw new Error("SC_SHINY_DATASET must be visium or slideseq");
const port = Number(process.env.SC_SHINY_PORT || 18766);
const debugPort = Number(process.env.SC_SHINY_CDP_PORT || 9232);
const chrome = process.env.CHROME_PATH ||
  "C:\\Program Files\\Google\\Chrome\\Application\\chrome.exe";
const tempRoot = require("os").tmpdir();
const userData = fs.mkdtempSync(path.join(tempRoot, "spatialcellchat-shiny-large-cdp-"));

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

async function waitForUrl(url, timeout = 60000) {
  const started = Date.now();
  while (Date.now() - started < timeout) {
    try {
      const response = await fetch(url);
      if (response.ok) return;
    } catch (_) {}
    await new Promise(resolve => setTimeout(resolve, 200));
  }
  throw new Error("Timed out waiting for " + url);
}

async function waitForExpression(client, expression, timeout = 180000) {
  const started = Date.now();
  while (Date.now() - started < timeout) {
    if (await evaluate(client, expression)) return;
    await new Promise(resolve => setTimeout(resolve, 300));
  }
  throw new Error("Timed out waiting for expression: " + expression);
}

async function click(client, selector) {
  const box = await evaluate(client, "(() => {" +
    "const node = document.querySelector(" + JSON.stringify(selector) + ");" +
    "if (!node) return null; const r = node.getBoundingClientRect();" +
    "return {x:r.left+r.width/2,y:r.top+r.height/2};})()");
  if (!box) throw new Error("Missing clickable element " + selector);
  await client.command("Input.dispatchMouseEvent", { type: "mousePressed", x: box.x, y: box.y, button: "left", clickCount: 1 });
  await client.command("Input.dispatchMouseEvent", { type: "mouseReleased", x: box.x, y: box.y, button: "left", clickCount: 1 });
}

async function dragSelect(client, x1, y1, x2, y2) {
  await client.command("Input.dispatchMouseEvent", { type: "mouseMoved", x: x1, y: y1 });
  await client.command("Input.dispatchMouseEvent", { type: "mousePressed", x: x1, y: y1, button: "left", clickCount: 1 });
  for (let i = 1; i <= 14; i += 1) {
    const t = i / 14;
    await client.command("Input.dispatchMouseEvent", { type: "mouseMoved", x: x1 + (x2 - x1) * t, y: y1 + (y2 - y1) * t, button: "left" });
    await new Promise(resolve => setTimeout(resolve, 20));
  }
  await client.command("Input.dispatchMouseEvent", { type: "mouseReleased", x: x2, y: y2, button: "left", clickCount: 1 });
}

async function selectedIds(client) {
  return evaluate(client, "(() => {" +
    "const text=document.querySelector('#selecting2')?.innerText||'';" +
    "return [...text.matchAll(/'([^']+)'/g)].map(value=>value[1]);})()");
}

async function targetPoints(client) {
  return evaluate(client, "(() => {" +
    "const gd=document.querySelector('#Celltype'); const traces=gd?gd._fullData||[]:[];" +
    "const points=[]; for (const trace of traces) {" +
    "if (!Array.isArray(trace.customdata)||!Array.isArray(trace.x)||!Array.isArray(trace.y)) continue;" +
    "for (let i=0;i<trace.customdata.length;i+=1) if (Number.isFinite(Number(trace.x[i]))&&Number.isFinite(Number(trace.y[i]))) points.push({id:String(trace.customdata[i]),x:Number(trace.x[i]),y:Number(trace.y[i])});}" +
    "if (!points.length) return null; points.sort((a,b)=>a.x-b.x||a.y-b.y||a.id.localeCompare(b.id));" +
    "return {first:points[0].id,second:points[points.length-1].id};})()");
}

async function pointPixel(client, cell) {
  return evaluate(client, "(() => {" +
    "const gd=document.querySelector('#Celltype'); if (!gd) return null;" +
    "const trace=(gd._fullData||[]).find(item=>Array.isArray(item.customdata)&&item.customdata.includes(" + JSON.stringify(cell) + "));" +
    "if (!trace) return null; const index=trace.customdata.indexOf(" + JSON.stringify(cell) + ");" +
    "const xa=gd._fullLayout.xaxis, ya=gd._fullLayout.yaxis, rect=gd.getBoundingClientRect();" +
    "return {x:rect.left+xa._offset+xa.l2p(trace.x[index]),y:rect.top+ya._offset+ya.l2p(trace.y[index])};})()");
}

async function main() {
  await waitForUrl("http://127.0.0.1:" + port + "/");
  const browser = spawn(chrome, ["--headless=new", "--disable-gpu", "--no-sandbox",
    "--remote-debugging-port=" + debugPort, "--user-data-dir=" + userData,
    "--window-size=1400,1000", "about:blank"], { windowsHide: true, stdio: "ignore" });
  let client;
  try {
    const list = await waitForUrl("http://127.0.0.1:" + debugPort + "/json/list")
      .then(() => fetch("http://127.0.0.1:" + debugPort + "/json/list").then(response => response.json()));
    const page = list.find(item => item.type === "page");
    if (!page) throw new Error("Chrome did not expose a page target");
    client = cdp(page.webSocketDebuggerUrl);
    await client.command("Page.enable");
    const started = Date.now();
    await client.command("Page.navigate", { url: "http://127.0.0.1:" + port + "/" });
    await waitForExpression(client, "!!document.querySelector('#Celltype .plotly')");
    await waitForExpression(client, "!!(document.querySelector('#Celltype')?._fullData||[]).length", 300000);
    await waitForExpression(client, "!!(document.querySelector('#LRExpr1')?._fullData||[]).length", 300000);
    await waitForExpression(client, "!!document.querySelector('#update')&&!!document.querySelector('#reset')");
    await waitForExpression(client, "((document.querySelector('#selecting2')?.innerText||'').match(/'[^']+'/g)||[]).length===" + expectedCells, 300000);
    const initial = await selectedIds(client);
    const expressionCells = await evaluate(client, "(() => {const gd=document.querySelector('#LRExpr1');return new Set((gd._fullData||[]).flatMap(trace=>Array.isArray(trace.customdata)?trace.customdata:[])).size;})()");
    if (expressionCells !== expectedCells) throw new Error("Expression plot covered " + expressionCells + " of " + expectedCells + " cells");
    const initialMs = Date.now() - started;
    const targets = await targetPoints(client);
    if (!targets) throw new Error("Could not identify target cells in Plotly data");
    const first = await pointPixel(client, targets.first);
    const second = await pointPixel(client, targets.second);
    if (!first || !second) throw new Error("Could not map target cells to browser pixels");

    const firstStarted = Date.now();
    await dragSelect(client, first.x - 10, first.y - 10, first.x + 10, first.y + 10);
    await waitForExpression(client, "(document.querySelector('#selecting')?.innerText||'').includes(" + JSON.stringify("'" + targets.first + "'") + ")");
    await click(client, "#update");
    await waitForExpression(client, "(() => {const text=document.querySelector('#selecting2')?.innerText||''; const n=(text.match(/'[^']+'/g)||[]).length; return n>0&&n<" + expectedCells + "&&text.includes(" + JSON.stringify("'" + targets.first + "'") + ");})()");
    const firstSelection = await selectedIds(client);
    const firstMs = Date.now() - firstStarted;

    await new Promise(resolve => setTimeout(resolve, 1000));
    const secondStarted = Date.now();
    await dragSelect(client, second.x - 10, second.y - 10, second.x + 10, second.y + 10);
    await waitForExpression(client, "(document.querySelector('#selecting')?.innerText||'').includes(" + JSON.stringify("'" + targets.second + "'") + ")");
    await click(client, "#update");
    await waitForExpression(client, "(() => {const text=document.querySelector('#selecting2')?.innerText||''; const n=(text.match(/'[^']+'/g)||[]).length; return n>0&&n<" + expectedCells + "&&text.includes(" + JSON.stringify("'" + targets.second + "'") + ")&&!text.includes(" + JSON.stringify("'" + targets.first + "',") + ");})()");
    const secondSelection = await selectedIds(client);
    const secondMs = Date.now() - secondStarted;

    await click(client, "#reset");
    await waitForExpression(client, "((document.querySelector('#selecting2')?.innerText||'').match(/'[^']+'/g)||[]).length===" + expectedCells, 300000);
    const reset = await selectedIds(client);
    const elapsed = Date.now() - started;
    if (initial.length !== expectedCells || !firstSelection.includes(targets.first) ||
        !secondSelection.includes(targets.second) || reset.length !== expectedCells ||
        firstSelection.some(cell => secondSelection.includes(cell)))
      throw new Error("unexpected selection sets initial=" + initial.length +
        " first=" + firstSelection.length + " second=" + secondSelection.length +
        " reset=" + reset.length + " targets=" + JSON.stringify(targets));
    const out = path.resolve(__dirname, "figures", "P4", "shiny-browser-summary-" + dataset + "-2026-11-10.csv");
    fs.writeFileSync(out, "dataset,cells,expression_cells,initial_count,first_count,second_count,reset_count,initial_ms,first_ms,second_ms,elapsed_ms,first_target,second_target,status\n" +
      dataset + "," + expectedCells + "," + expressionCells + "," + initial.length + "," + firstSelection.length + "," + secondSelection.length + "," + reset.length + "," + initialMs + "," + firstMs + "," + secondMs + "," + elapsed + "," + JSON.stringify(targets.first) + "," + JSON.stringify(targets.second) + ",PASS\n");
    console.log("[PASS] " + dataset + " Shiny browser cells=" + initial.length + " first=" + firstSelection.length + " second=" + secondSelection.length + " reset=" + reset.length + " elapsed=" + elapsed + "ms");
  } finally {
    if (client) client.close();
    if (browser.exitCode === null) { browser.kill(); await once(browser, "exit"); }
    const resolvedTempRoot = path.resolve(tempRoot);
    const resolvedUserData = path.resolve(userData);
    if (path.dirname(resolvedUserData) !== resolvedTempRoot || !path.basename(resolvedUserData).startsWith("spatialcellchat-shiny-large-cdp-"))
      throw new Error("Refusing to remove an unexpected Chrome profile path");
    fs.rmSync(resolvedUserData, { recursive: true, force: true });
  }
}

main().catch(error => { console.error(error.stack || error); process.exitCode = 1; });
