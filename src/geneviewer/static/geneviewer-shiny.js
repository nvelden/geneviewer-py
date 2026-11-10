// The Shiny for Python output binding of geneviewer: draws each chart's spec
// with geneviewer.newChart() (geneviewer.min.js) and sends its events to
// Shiny inputs. Written after the R package's htmlwidgets binding
// (inst/htmlwidgets/geneviewerwidget.js), so both send the same inputs.
// Not written by scripts/update_js.py.
(function () {
  // The events sent to Shiny, as the inputs <output id>_<name>
  var shinyEvents = ["click", "dblclick", "select", "selection", "drag", "moves", "zoom", "rendered"];
  // The chart of each output element. A chart lives as long as its element,
  // so removing a Shiny output frees it.
  var charts = new WeakMap();

  var binding = new Shiny.OutputBinding();

  binding.find = function (scope) {
    return $(scope).find(".geneviewer-output");
  };

  // value is { spec }, or null when the render function returned None
  binding.renderValue = function (el, value) {
    var chart = charts.get(el);
    if (!value) {
      if (chart) chart.destroy();
      charts.delete(el);
      return;
    }
    // A redraw from Shiny gives the whole chart again, so it replaces the
    // chart's spec (setSpec), which update() would merge into
    if (chart) {
      chart.setSpec(value.spec);
      return;
    }
    var on = {};
    shinyEvents.forEach(function (name) {
      on[name] = function (event) {
        // Priority "event": the same event twice still reaches the server
        Shiny.setInputValue(el.id + "_" + name, event.detail, { priority: "event" });
      };
    });
    // geneviewer.css comes with the output's HTML dependency
    charts.set(el, geneviewer.newChart(el, value.spec, { on: on, server: true, injectStyles: false }));
  };

  binding.resize = function (el) {
    var chart = charts.get(el);
    if (chart) chart.resize();
  };

  Shiny.outputBindings.register(binding, "geneviewer.output");

  // zoom_to(): the chart of the output, found when the message arrives
  Shiny.addCustomMessageHandler("geneviewer-zoom", function (message) {
    var el = document.getElementById(message.id);
    var chart = el && charts.get(el);
    if (chart) chart.zoomTo(message.range);
  });
})();
