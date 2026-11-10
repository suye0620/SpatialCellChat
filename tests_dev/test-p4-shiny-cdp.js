/* Browser interaction smoke for the Shiny LR selector. */
const fs = require("fs");
const path = require("path");
const { spawn } = require("child_process");
const { once } = require("events");

const port = Number(process.env.SC_SHINY_PORT || 18765);
const debugPort = Number(process.env.SC_SHINY_CDP_PORT || 9230);
const chrome = process.env.CHROME_PATH ||
  "C:\\Program Files\\Google\\Chrome\\Application\\chrome.exe";
const tempRoot = require("os").tmpdir();
const userData = fs.mkdtempSync(path.join(tempRoot, "spatialcellchat-shiny-cdp-"));

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
    if (message.id && pending.has(message.id)) {
      const waiter = pending.get(message.id); pending.delete(message.id);
      if (message.error) waiter.reject(new Error(message.error.message));
      else waiter.resolve(message.result);
    }
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

async function waitForUrl(url, timeout = 30000) {
  const started = Date.now();
  while (Date.now() - started < timeout) {
    try {
      const response = await fetch(url);
      if (response.ok) return;
    } catch (_) {}
    await new Promise(resolve => setTimeout(resolve, 100));
  }
  throw new Error(`Timed out waiting for ${url}`);
}

async function evaluate(client, expression) {
  const result = await client.command("Runtime.evaluate", {
    expression, awaitPromise: true, returnByValue: true
  });
  if (result.exceptionDetails)
    throw new Error(result.exceptionDetails.text || "browser evaluation failed");
  return result.result.value;
}

async function waitForExpression(client, expression, timeout = 30000) {
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
  await client.command("Input.dispatchMouseEvent", {
    type: "mousePressed", x: box.x, y: box.y, button: "left", clickCount: 1
  });
  await client.command("Input.dispatchMouseEvent", {
    type: "mouseReleased", x: box.x, y: box.y, button: "left", clickCount: 1
  });
}

async function dragSelect(client, x1, y1, x2, y2) {
  await client.command("Input.dispatchMouseEvent", {
    type: "mouseMoved", x: x1, y: y1
  });
  await client.command("Input.dispatchMouseEvent", {
    type: "mousePressed", x: x1, y: y1, button: "left", clickCount: 1
  });
  for (let i = 1; i <= 12; i++) {
    const t = i / 12;
    await client.command("Input.dispatchMouseEvent", {
      type: "mouseMoved", x: x1 + (x2 - x1) * t, y: y1 + (y2 - y1) * t,
      button: "left"
    });
    await new Promise(resolve => setTimeout(resolve, 20));
  }
  await client.command("Input.dispatchMouseEvent", {
    type: "mouseReleased", x: x2, y: y2, button: "left", clickCount: 1
  });
}

async function selectionIds(client) {
  return evaluate(client, `(() => {
    const text = document.querySelector('#selecting2')?.innerText || '';
    return (text.match(/'c\\d+'/g) || []).map(value => value.slice(1, -1));
  })()`);
}

async function plotBox(client) {
  return evaluate(client, `(() => {
    const node = document.querySelector('#Celltype .plotly');
    if (!node) return null;
    const r = node.getBoundingClientRect();
    return {left: r.left, top: r.top, width: r.width, height: r.height};
  })()`);
}

async function pointPixel(client, cell) {
  return evaluate(client, `(() => {
    const gd = document.querySelector('#Celltype');
    const trace = (gd?._fullData || []).find(item =>
      Array.isArray(item.customdata) && item.customdata.includes(${JSON.stringify(cell)}));
    if (!gd) return {error: 'plot missing'};
    if (!trace) return {error: 'cell missing', traces: (gd._fullData || []).map(item => ({
      type: item.type, customdata: item.customdata && item.customdata.slice(0, 5),
      text: item.text && item.text.slice(0, 2)
    }))};
    const index = trace.customdata.indexOf(${JSON.stringify(cell)});
    const xa = gd._fullLayout.xaxis;
    const ya = gd._fullLayout.yaxis;
    const rect = gd.querySelector('.plotly').getBoundingClientRect();
    return {
      x: rect.left + xa._offset + xa.l2p(trace.x[index]),
      y: rect.top + ya._offset + ya.l2p(trace.y[index])
    };
  })()`);
}

async function main() {
  await waitForUrl(`http://127.0.0.1:${port}/`);
  const browser = spawn(chrome, [
    "--headless=new", "--disable-gpu", "--no-sandbox",
    `--remote-debugging-port=${debugPort}`, `--user-data-dir=${userData}`,
    "--window-size=1400,1000", "about:blank"
  ], { windowsHide: true, stdio: "ignore" });
  try {
    const pages = await waitForUrl(`http://127.0.0.1:${debugPort}/json/list`)
      .then(() => fetch(`http://127.0.0.1:${debugPort}/json/list`).then(r => r.json()));
    const page = pages.find(item => item.type === "page");
    if (!page) throw new Error("Chrome did not expose a page target");
    const client = cdp(page.webSocketDebuggerUrl);
    await client.command("Page.enable");
    await client.command("Page.navigate", { url: `http://127.0.0.1:${port}/` });
    await waitForExpression(client, `!!document.querySelector('#Celltype .plotly')`);
    try {
      await waitForExpression(client, `(() => {
        const gd = document.querySelector('#Celltype');
        return !!(gd && gd._fullData && gd._fullData.length);
      })()`);
    } catch (error) {
      const debug = await evaluate(client, `(() => {
        const gd = document.querySelector('#Celltype');
        return {plotly: typeof window.Plotly, classes: Array.from(document.querySelectorAll('#Celltype *'))
            .map(node => String(node.className || '')).filter(name => name.includes('plotly')).slice(0, 10),
          keys: gd ? Object.keys(gd).slice(0, 20) : [],
          fullData: gd && gd._fullData ? gd._fullData.length : null,
          html: document.querySelector('#Celltype')?.innerHTML.slice(0, 500)};
      })()`);
      throw new Error(`${error.message}; debug=${JSON.stringify(debug)}`);
    }
    await waitForExpression(client, `!!document.querySelector('#update') && !!document.querySelector('#reset')`);
    await waitForExpression(client, `(() => {
      const text = document.querySelector('#selecting2')?.innerText || '';
      return (text.match(/'c\\d+'/g) || []).length === 20;
    })()`);
    const initial = await selectionIds(client);
    const box = await plotBox(client);
    const firstPoint = await pointPixel(client, "c1");
    if (!box || !firstPoint || firstPoint.error)
      throw new Error(`Plotly spatial plot has no usable browser geometry: ${JSON.stringify({box, firstPoint})}`);
    await dragSelect(client, firstPoint.x - 12, firstPoint.y - 12,
                     firstPoint.x + 12, firstPoint.y + 12);
    await waitForExpression(client, `document.querySelector('#selecting')?.innerText.includes("'c1'")`);
    await click(client, "#update");
    await waitForExpression(client, `(() => {
      const text = document.querySelector('#selecting2')?.innerText || '';
      const n = (text.match(/'c\\d+'/g) || []).length;
      return n > 0 && n < 20;
    })()`);
    const first = await selectionIds(client);
    await new Promise(resolve => setTimeout(resolve, 1500));
    const secondPoint = await pointPixel(client, "c13");
    if (!secondPoint) throw new Error("Plotly spatial plot lost the second target cell");
    await dragSelect(client, secondPoint.x - 12, secondPoint.y - 12,
                     secondPoint.x + 12, secondPoint.y + 12);
    await waitForExpression(client, `document.querySelector('#selecting')?.innerText.includes("'c13'")`);
    await click(client, "#update");
    await waitForExpression(client, `(() => {
      const text = document.querySelector('#selecting2')?.innerText || '';
      return text.includes("'c13'") && !text.includes("'c1',");
    })()`);
    const second = await selectionIds(client);
    await click(client, "#reset");
    await waitForExpression(client, `(() => {
      const text = document.querySelector('#selecting2')?.innerText || '';
      return (text.match(/'c\\d+'/g) || []).length === 20;
    })()`);
    const reset = await selectionIds(client);
    const disjoint = first.length > 0 && second.length > 0 &&
      first.every(cell => !second.includes(cell));
    if (!(initial.length === 20 && first.length > 0 && first.length < 20 &&
          second.length > 0 && second.length < 20 && disjoint && reset.length === 20))
      throw new Error(`unexpected selection sets initial=${initial} first=${first} second=${second} reset=${reset}`);
    const summary = `dataset,initial_count,first_selection,second_selection,reset_count,status\n` +
      `simulated,${initial.length},${first.join("|")},${second.join("|")},${reset.length},PASS\n`;
    fs.writeFileSync(path.resolve(__dirname, "figures", "P4",
      "shiny-browser-summary-2026-11-10.csv"), summary);
    console.log(`[PASS] shiny browser selection initial=${initial.length} first=${first.join("|")} second=${second.join("|")} reset=${reset.length}`);
    client.close();
  } finally {
    if (browser.exitCode === null) {
      browser.kill();
      await once(browser, "exit");
    }
    const resolvedTempRoot = path.resolve(tempRoot);
    const resolvedUserData = path.resolve(userData);
    if (path.dirname(resolvedUserData) !== resolvedTempRoot ||
        !path.basename(resolvedUserData).startsWith("spatialcellchat-shiny-cdp-"))
      throw new Error("Refusing to remove an unexpected Chrome profile path");
    fs.rmSync(resolvedUserData, { recursive: true, force: true });
  }
}

main().catch(error => { console.error(error.stack || error); process.exitCode = 1; });
