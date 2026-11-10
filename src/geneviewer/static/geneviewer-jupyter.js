// The Jupyter widget of geneviewer (anywidget): draws the chart's spec with
// newChart() and sends its events to Python. geneviewer.jupyter puts this
// file after the library (geneviewer.standalone.mjs) in one module and gives
// it the library's newChart. Written after the Shiny binding
// (geneviewer-shiny.js), so both report the same events.
// Not written by scripts/update_js.py.
const __geneviewerRender = (function (newChart) {
  // The events sent to Python, as the widget's attributes of the same name
  const events = ["click", "dblclick", "select", "selection", "drag", "moves", "zoom", "rendered"];
  const cssSize = (value, fallback) => {
    if (value === undefined || value === null) return fallback;
    return typeof value === "number" ? value + "px" : String(value);
  };

  return function render({ model, el }) {
    // The chart fills its element, which fills a frame of the spec's size
    const frame = document.createElement("div");
    frame.className = "geneviewer-jupyter";
    const holder = document.createElement("div");
    // The chart is named after its element. Each widget loads this module
    // anew, so the library's own count of charts starts at one in each:
    // a name of its own keeps the charts of a notebook apart
    holder.id = "geneviewer-widget-" + Math.random().toString(36).slice(2, 12);
    frame.appendChild(holder);
    el.appendChild(frame);
    const size = () => {
      const spec = model.get("spec") || {};
      frame.style.width = cssSize(spec.width, "100%");
      frame.style.height = cssSize(spec.height, "400px");
    };
    size();

    const on = {};
    events.forEach((name) => {
      on[name] = (event) => {
        // As a message, so that the same event twice reaches Python twice.
        // Through JSON: the detail is plain data, as Shiny gets it
        const detail = event.detail === undefined ? null : JSON.parse(JSON.stringify(event.detail));
        model.send({ type: "event", name, detail });
      };
    });
    const chart = newChart(holder, model.get("spec"), { on });

    // A new chart from Python replaces the spec (update() would merge into it)
    const redraw = () => {
      size();
      chart.setSpec(model.get("spec"));
    };
    model.on("change:spec", redraw);
    // zoom_to()
    const message = (content) => {
      if (content && content.type === "zoom") chart.zoomTo(content.range);
    };
    model.on("msg:custom", message);

    let width = holder.clientWidth;
    const sizes = new ResizeObserver(() => {
      if (holder.clientWidth === width) return;
      width = holder.clientWidth;
      chart.resize();
    });
    sizes.observe(holder);

    return () => {
      sizes.disconnect();
      model.off("change:spec", redraw);
      model.off("msg:custom", message);
      chart.destroy();
      frame.remove();
    };
  };
})(__GENEVIEWER_NEW_CHART__);

export default { render: __geneviewerRender };
