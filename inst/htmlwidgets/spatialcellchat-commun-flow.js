(function () {
  "use strict";

  function finite(value, fallback) {
    return Number.isFinite(value) ? value : fallback;
  }

  function clamp(value, lower, upper) {
    return Math.max(lower, Math.min(upper, value));
  }

  HTMLWidgets.widget({
    name: "spatialcellchat-commun-flow",
    type: "output",
    factory: function (el, width, height) {
      var state = {
        root: null,
        canvas: null,
        trailCanvas: null,
        staticCanvas: null,
        ctx: null,
        trailCtx: null,
        staticCtx: null,
        data: null,
        particles: [],
        running: false,
        rafId: null,
        lastTime: null,
        frameSeconds: 0,
        rngState: 1,
        scale: null,
        viewBounds: null,
        fieldBounds: null,
        maxMagnitude: 0,
        controls: null,
        tooltip: null,
        status: null,
        destroyed: false
      };

      function random() {
        state.rngState = (Math.imul(1664525, state.rngState) + 1013904223) >>> 0;
        return state.rngState / 4294967296;
      }

      function field() {
        return state.data.field;
      }

      function fieldBounds() {
        if (state.fieldBounds) return state.fieldBounds;
        var grid = field();
        state.fieldBounds = {
          xmin: Math.min.apply(null, grid.x),
          xmax: Math.max.apply(null, grid.x),
          ymin: Math.min.apply(null, grid.y),
          ymax: Math.max.apply(null, grid.y)
        };
        return state.fieldBounds;
      }

      function index(i, j) {
        return j * field().x.length + i;
      }

      function bracket(values, value) {
        var n = values.length;
        if (!n) return null;
        if (n === 1) return { i: 0, next: 0, t: 0 };
        var ascending = values[n - 1] >= values[0];
        var v = value;
        if (ascending) v = clamp(v, values[0], values[n - 1]);
        else v = clamp(v, values[n - 1], values[0]);
        var lo = 0;
        var hi = n - 1;
        while (hi - lo > 1) {
          var mid = Math.floor((lo + hi) / 2);
          if (ascending ? values[mid] <= v : values[mid] >= v) lo = mid;
          else hi = mid;
        }
        var denominator = values[lo === hi ? lo : hi] - values[lo];
        var t = denominator === 0 ? 0 : (v - values[lo]) / denominator;
        return { i: lo, next: hi, t: clamp(t, 0, 1) };
      }

      function sample(x, y) {
        var grid = field();
        var bx = bracket(grid.x, x);
        var by = bracket(grid.y, y);
        if (!bx || !by) return null;
        var i00 = index(bx.i, by.i);
        var i10 = index(bx.next, by.i);
        var i01 = index(bx.i, by.next);
        var i11 = index(bx.next, by.next);
        if (!grid.valid[i00] || !grid.valid[i10] ||
            !grid.valid[i01] || !grid.valid[i11]) return null;
        var tx = bx.t;
        var ty = by.t;
        var u0 = grid.u[i00] * (1 - tx) + grid.u[i10] * tx;
        var u1 = grid.u[i01] * (1 - tx) + grid.u[i11] * tx;
        var v0 = grid.v[i00] * (1 - tx) + grid.v[i10] * tx;
        var v1 = grid.v[i01] * (1 - tx) + grid.v[i11] * tx;
        return {
          u: u0 * (1 - ty) + u1 * ty,
          v: v0 * (1 - ty) + v1 * ty
        };
      }


      function pixel(x, y) {
        var bounds = state.viewBounds || state.data.bounds;
        var dx = bounds.xmax - bounds.xmin || 1;
        var dy = bounds.ymax - bounds.ymin || 1;
        var fit = Math.min(state.canvas.width / dx, state.canvas.height / dy);
        var offsetX = (state.canvas.width - dx * fit) / 2;
        var offsetY = (state.canvas.height - dy * fit) / 2;
        return {
          x: offsetX + (x - bounds.xmin) * fit,
          y: state.canvas.height - offsetY - (y - bounds.ymin) * fit
        };
      }

      function hasCurrent() {
        var grid = field();
        for (var i = 0; i < grid.valid.length; i += 1) {
          if (grid.valid[i] && (Math.abs(grid.u[i]) > 0 || Math.abs(grid.v[i]) > 0)) return true;
        }
        return false;
      }

      function calculateMaxMagnitude() {
        var grid = field();
        var maximum = 0;
        for (var i = 0; i < grid.u.length; i += 1) {
          if (grid.valid[i]) maximum = Math.max(maximum, Math.hypot(grid.u[i], grid.v[i]));
        }
        var cells = state.data.cell;
        for (var c = 0; c < cells.magnitude.length; c += 1) {
          maximum = Math.max(maximum, finite(cells.magnitude[c], 0));
        }
        return maximum;
      }

      function seedParticle() {
        var grid = field();
        var bounds = fieldBounds();
        for (var attempt = 0; attempt < 40; attempt += 1) {
          var x = bounds.xmin + random() * (bounds.xmax - bounds.xmin || 1);
          var y = bounds.ymin + random() * (bounds.ymax - bounds.ymin || 1);
          var velocity = sample(x, y);
          if (velocity && (Math.abs(velocity.u) > 0 || Math.abs(velocity.v) > 0)) {
            return { x: x, y: y, px: x, py: y, age: random() * (state.data.options.trailLength || 0) };
          }
        }
        for (var j = 0; j < grid.y.length; j += 1) {
          for (var i = 0; i < grid.x.length; i += 1) {
            var at = index(i, j);
            if (grid.valid[at] && (Math.abs(grid.u[at]) > 0 || Math.abs(grid.v[at]) > 0)) {
              return { x: grid.x[i], y: grid.y[j], px: grid.x[i], py: grid.y[j], age: 0 };
            }
          }
        }
        return null;
      }

      function resetParticles() {
        state.particles = [];
        if (!hasCurrent()) return;
        var count = state.data.options.particleCount;
        for (var i = 0; i < count; i += 1) {
          var particle = seedParticle();
          if (particle) state.particles.push(particle);
        }
      }



      function drawCells() {
        var ctx = state.staticCtx;
        var cells = state.data.cell;
        var style = state.data.style || {};
        var dpr = state.scale ? state.scale.pixelRatio : (window.devicePixelRatio || 1);
        ctx.clearRect(0, 0, state.staticCanvas.width, state.staticCanvas.height);
        state.maxMagnitude = calculateMaxMagnitude();
        var radius = Math.max(3, finite(style["point.size"], 1.4) * 2.2) * dpr;
        var fillAlpha = clamp(finite(style["image.alpha"], 0.32), 0.05, 1);
        for (var i = 0; i < cells.x.length; i += 1) {
          var point = pixel(cells.x[i], cells.y[i]);
          ctx.save();
          ctx.globalAlpha = fillAlpha;
          ctx.fillStyle = cells.color[i] || "#54748a";
          ctx.beginPath();
          ctx.arc(point.x, point.y, radius, 0, Math.PI * 2);
          ctx.fill();
          ctx.globalAlpha = 1;
          ctx.strokeStyle = "#ffffff";
          ctx.lineWidth = Math.max(1.1, 1.3 * dpr);
          ctx.stroke();
          ctx.restore();
        }
      }

      function drawParticles() {
        var ctx = state.trailCtx;
        var options = state.data.options;
        var trail = Math.max(0, options.trailLength || 0);
        var fade = trail === 0 ? 1 : clamp(0.32 - trail * 0.014, 0.04, 0.32);
        ctx.save();
        if (fade >= 1) ctx.clearRect(0, 0, state.trailCanvas.width, state.trailCanvas.height);
        else {
          ctx.fillStyle = "rgba(247, 249, 246, " + fade + ")";
          ctx.fillRect(0, 0, state.trailCanvas.width, state.trailCanvas.height);
        }
        for (var p = 0; p < state.particles.length; p += 1) {
          var particle = state.particles[p];
          var velocity = sample(particle.x, particle.y);
          if (!velocity || (!velocity.u && !velocity.v) || particle.age > trail + 1) {
            var replacement = seedParticle();
            if (replacement) state.particles[p] = replacement;
            continue;
          }
          var magnitude = Math.hypot(velocity.u, velocity.v);
          var brightness = state.maxMagnitude ? Math.sqrt(magnitude / state.maxMagnitude) : 0;
          var alpha = 0.24 + brightness * 0.55;
          var from = pixel(particle.x, particle.y);
          var maxStep = state.scale.step * 0.36;
          var factor = options.particleSpeed * state.frameSeconds * maxStep /
            (state.maxMagnitude || 1);
          particle.px = particle.x;
          particle.py = particle.y;
          particle.x += velocity.u * factor;
          particle.y += velocity.v * factor;
          particle.age += state.frameSeconds * options.particleSpeed;
          var to = pixel(particle.x, particle.y);
          if (sample(particle.x, particle.y) === null ||
              particle.x < fieldBounds().xmin || particle.x > fieldBounds().xmax ||
              particle.y < fieldBounds().ymin || particle.y > fieldBounds().ymax) {
            var restarted = seedParticle();
            if (restarted) state.particles[p] = restarted;
            continue;
          }
          ctx.lineWidth = (0.8 + brightness * 1.1) * state.scale.pixelRatio;
          ctx.lineCap = "round";
          ctx.beginPath();
          ctx.moveTo(from.x, from.y);
          ctx.lineTo(to.x, to.y);
          ctx.stroke();
          ctx.fillStyle = "rgba(19, 117, 145, " + Math.min(0.9, alpha + 0.18).toFixed(3) + ")";
          ctx.beginPath();
          ctx.arc(to.x, to.y, (1.5 + brightness * 1.1) * state.scale.pixelRatio, 0, Math.PI * 2);
          ctx.fill();
        }
        ctx.restore();
      }

      function frame(timestamp) {
        if (!state.running || state.destroyed) return;
        var previous = state.lastTime === null ? timestamp : state.lastTime;
        state.frameSeconds = clamp((timestamp - previous) / 1000, 0, 0.05);
        state.lastTime = timestamp;
        drawParticles();
        state.rafId = window.requestAnimationFrame(frame);
      }

      function setRunning(value) {
        state.running = value;
        state.lastTime = null;
        if (!value && state.rafId !== null) {
          window.cancelAnimationFrame(state.rafId);
          state.rafId = null;
        }
        if (state.controls && state.controls.button) {
          state.controls.button.textContent = value ? "Pause" : "Play";
          state.controls.button.setAttribute("aria-label", value ? "Pause animation" : "Play animation");
        }
        if (value && state.rafId === null && hasCurrent()) {
          state.rafId = window.requestAnimationFrame(function (timestamp) {
            state.rafId = null;
            frame(timestamp);
          });
        }
      }

      function buildLegend() {
        var legend = document.createElement("div");
        legend.className = "sc-flow-legend";
        var heading = document.createElement("span");
        heading.className = "sc-flow-legend-title";
        var metadata = state.data.metadata || {};
        heading.textContent = metadata.legendTitle ||
          (metadata.pattern === "incoming" ? "Targets" : "Sources");
        legend.appendChild(heading);
        var colors = state.data.style && state.data.style.colors ? state.data.style.colors : {};
        var labels = state.data.cell.label.filter(function (value, i, array) {
          return array.indexOf(value) === i;
        });
        labels.forEach(function (label) {
          var item = document.createElement("span");
          item.className = "sc-flow-legend-item";
          var swatch = document.createElement("span");
          swatch.className = "sc-flow-legend-swatch";
          swatch.style.backgroundColor = colors[label] || "#718096";
          var text = document.createElement("span");
          text.textContent = label;
          item.appendChild(swatch);
          item.appendChild(text);
          legend.appendChild(item);
        });
        var direction = document.createElement("span");
        direction.className = "sc-flow-legend-item sc-flow-legend-direction";
        direction.textContent = "particles = signal-flow direction";
        legend.appendChild(direction);
        return legend;
      }

      function buildRoot() {
        var root = document.createElement("div");
        root.className = "sc-flow-root";
        var header = document.createElement("div");
        header.className = "sc-flow-header";
        var title = document.createElement("div");
        title.className = "sc-flow-title";
        title.textContent = state.data.options.title || "Communication flow";
        var semantics = document.createElement("div");
        semantics.className = "sc-flow-semantics";
        var metadata = state.data.metadata || {};
        semantics.textContent = (metadata.pattern === "incoming" ? "Incoming" : "Outgoing") +
          " · particle motion = signal flow";
        header.appendChild(title);
        header.appendChild(semantics);
        root.appendChild(header);

        var controls = document.createElement("div");
        controls.className = "sc-flow-controls";
        var button = document.createElement("button");
        button.type = "button";
        button.className = "sc-flow-button";
        button.addEventListener("click", function () { setRunning(!state.running); });
        controls.appendChild(button);

        function addRange(labelText, min, max, step, value, callback) {
          var label = document.createElement("label");
          label.className = "sc-flow-control";
          var text = document.createElement("span");
          text.textContent = labelText;
          var input = document.createElement("input");
          input.type = "range";
          input.min = min;
          input.max = max;
          input.step = step;
          input.value = value;
          var output = document.createElement("output");
          output.className = "sc-flow-control-value";
          output.textContent = value;
          input.addEventListener("input", function () {
            output.textContent = input.value;
            callback(Number(input.value));
          });
          label.appendChild(text);
          label.appendChild(input);
          label.appendChild(output);
          controls.appendChild(label);
          return input;
        }

        var speed = addRange("Visual speed", 0.05, 3, 0.05,
          state.data.options.particleSpeed, function (value) {
            state.data.options.particleSpeed = value;
          });
        var count = addRange("Particles", 1, 5000, 1,
          state.data.options.particleCount, function (value) {
            state.data.options.particleCount = value;
            resetParticles();
            updateStatus();
          });
        var opacity = addRange("Point opacity", 0.05, 1, 0.05,
          state.data.style["image.alpha"], function (value) {
            state.data.style["image.alpha"] = value;
            drawCells();
          });
        var zoom = addRange("Zoom", 1, 4, 0.1,
          state.data.options.zoom || 1, function (value) {
            state.data.options.zoom = value;
            updateViewBounds();
            clearTrails();
            drawCells();
          });
        state.controls = {
          button: button, speed: speed, count: count, opacity: opacity, zoom: zoom
        };
        root.appendChild(controls);

        var plot = document.createElement("div");
        plot.className = "sc-flow-plot";
        plot.setAttribute("role", "img");
        plot.setAttribute("aria-label", state.data.options.title || "Communication flow plot");
        state.staticCanvas = document.createElement("canvas");
        state.trailCanvas = document.createElement("canvas");
        state.canvas = state.trailCanvas;
        state.staticCanvas.className = "sc-flow-canvas sc-flow-cells";
        state.trailCanvas.className = "sc-flow-canvas sc-flow-trails";
        plot.appendChild(state.staticCanvas);
        plot.appendChild(state.trailCanvas);
        state.trailCanvas.addEventListener("wheel", function (event) {
          event.preventDefault();
          var currentZoom = finite(state.data.options.zoom, 1);
          var nextZoom = clamp(currentZoom * Math.exp(-event.deltaY * 0.0015), 1, 4);
          state.data.options.zoom = Math.round(nextZoom * 10) / 10;
          if (state.controls && state.controls.zoom) {
            state.controls.zoom.value = state.data.options.zoom;
            state.controls.zoom.dispatchEvent(new Event("input"));
          }
        }, { passive: false });
        state.tooltip = document.createElement("div");
        state.tooltip.className = "sc-flow-tooltip";
        state.tooltip.hidden = true;
        plot.appendChild(state.tooltip);
        root.appendChild(plot);
        root.appendChild(buildLegend());
        state.status = document.createElement("div");
        state.status.className = "sc-flow-status";
        root.appendChild(state.status);

        state.trailCanvas.addEventListener("pointermove", function (event) {
          var rect = state.trailCanvas.getBoundingClientRect();
          var px = (event.clientX - rect.left) * state.trailCanvas.width / rect.width;
          var py = (event.clientY - rect.top) * state.trailCanvas.height / rect.height;
          var nearest = -1;
          var nearestDistance = Infinity;
          var cells = state.data.cell;
          for (var i = 0; i < cells.x.length; i += 1) {
            var point = pixel(cells.x[i], cells.y[i]);
            var distance = Math.hypot(point.x - px, point.y - py);
            if (distance < nearestDistance) {
              nearest = i;
              nearestDistance = distance;
            }
          }
          if (nearest >= 0 && nearestDistance <= 18 * state.scale.pixelRatio) {
            var cell = state.data.cell;
            var mass = Number.isFinite(cell.mass[nearest]) ?
              " · mass " + cell.mass[nearest].toPrecision(4) : "";
            state.tooltip.textContent = cell.id[nearest] + " · " + cell.label[nearest] +
              " · current " + cell.magnitude[nearest].toPrecision(4) + mass;
            state.tooltip.style.left = Math.min(px / state.scale.pixelRatio + 12,
              rect.width - 230) + "px";
            state.tooltip.style.top = Math.max(py / state.scale.pixelRatio - 30, 6) + "px";
            state.tooltip.hidden = false;
          } else {
            state.tooltip.hidden = true;
          }
        });
        state.trailCanvas.addEventListener("pointerleave", function () {
          state.tooltip.hidden = true;
        });
        return root;
      }

      function updateStatus() {
        if (!state.status) return;
        state.status.textContent = hasCurrent() ?
          state.data.cell.id.length + " cells · " + state.particles.length + " particles" :
          "No non-zero net current in this field";
      }

      function clearTrails() {
        if (state.trailCtx && state.trailCanvas)
          state.trailCtx.clearRect(0, 0, state.trailCanvas.width, state.trailCanvas.height);
      }

      function updateViewBounds() {
        if (!state.scale || !state.data) return;
        var bounds = state.data.bounds;
        var xRange = bounds.xmax - bounds.xmin || 1;
        var yRange = bounds.ymax - bounds.ymin || 1;
        var padding = Math.max(Math.max(xRange, yRange) * 0.14, state.scale.step * 0.8);
        var zoom = clamp(finite(state.data.options.zoom, 1), 1, 4);
        var centerX = (bounds.xmin + bounds.xmax) / 2;
        var centerY = (bounds.ymin + bounds.ymax) / 2;
        var halfX = (xRange + 2 * padding) / (2 * zoom);
        var halfY = (yRange + 2 * padding) / (2 * zoom);
        state.viewBounds = {
          xmin: centerX - halfX,
          xmax: centerX + halfX,
          ymin: centerY - halfY,
          ymax: centerY + halfY
        };
      }

      function resizeCanvases(newWidth, newHeight) {
        if (!state.root) return;
        var plot = state.canvas.parentElement;
        var plotWidth = Math.max(240, plot.clientWidth || newWidth || 640);
        var bounds = state.data.bounds;
        var ratio = (bounds.ymax - bounds.ymin) / (bounds.xmax - bounds.xmin || 1);
        var plotHeight = newHeight && newHeight > 80 ? newHeight : clamp(plotWidth * ratio, 280, 680);
        plot.style.height = plotHeight + "px";
        var dpr = window.devicePixelRatio || 1;
        state.scale = { step: Math.min(
          Math.abs(field().x[1] - field().x[0] || 1),
          Math.abs(field().y[1] - field().y[0] || 1)
        ), pixelRatio: dpr };
        state.fieldBounds = null;
        updateViewBounds();
        [state.trailCanvas, state.staticCanvas].forEach(function (canvas) {
          canvas.style.width = plotWidth + "px";
          canvas.style.height = plotHeight + "px";
          canvas.width = Math.round(plotWidth * dpr);
          canvas.height = Math.round(plotHeight * dpr);
        });
        state.ctx = state.trailCanvas.getContext("2d");
        state.trailCtx = state.ctx;
        state.staticCtx = state.staticCanvas.getContext("2d");
        drawCells();
        updateStatus();
      }

      function renderValue(x) {
        destroy();
        state.destroyed = false;
        state.data = x;
        state.rngState = (Number(x.options.seed) >>> 0) || 1;
        state.root = buildRoot();
        el.replaceChildren(state.root);
        resetParticles();
        resizeCanvases(width, height);
        setRunning(hasCurrent());
      }

      function resize(newWidth, newHeight) {
        resizeCanvases(newWidth, newHeight);
      }

      function destroy() {
        state.destroyed = true;
        state.running = false;
        if (state.rafId !== null) {
          window.cancelAnimationFrame(state.rafId);
          state.rafId = null;
        }
        state.lastTime = null;
        state.root = null;
        state.controls = null;
        state.tooltip = null;
        state.status = null;
        state.fieldBounds = null;
        if (el) el.replaceChildren();
      }

      el.__spatialCellChatFlow = {
        pause: function () { setRunning(false); },
        resume: function () { setRunning(true); },
        isRunning: function () { return state.running; },
        particleSnapshot: function () {
          return state.particles.map(function (particle) { return [particle.x, particle.y]; });
        },
        particleCount: function () { return state.particles.length; }
      };

      return {
        renderValue: renderValue,
        resize: resize,
        destroy: destroy
      };
    }
  });
}());
