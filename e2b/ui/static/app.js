// Two behaviours only. No framework, no build step, no key material.
(function () {
  "use strict";

  // 1. Suppress double submission. The server also de-duplicates by one-shot token;
  //    this only removes the visual invitation to click twice.
  document.querySelectorAll("form[data-once]").forEach(function (form) {
    form.addEventListener("submit", function () {
      var btn = form.querySelector("button[type=submit]");
      if (btn) {
        if (btn.dataset.submitted === "1") { return false; }
        btn.dataset.submitted = "1";
        btn.disabled = true;
        btn.textContent = "Submitting…";
      }
    });
  });

  // 2. Preset prefill. Fills the form fields and nothing else; the fields stay editable.
  document.querySelectorAll("[data-preset]").forEach(function (b) {
    b.addEventListener("click", function () {
      var set = function (id, val) { var el = document.getElementById(id); if (el) el.value = val || ""; };
      set("raw_indication", b.dataset.raw);
      set("desired_effect", b.dataset.effect);
      set("moa", b.dataset.moa);
      set("modality", b.dataset.modality);
    });
  });

  // 3. Poll the real job fragment while the job is non-terminal. Replaces the fragment
  //    with whatever the server says. Never animates, never interpolates a percentage.
  var node = document.getElementById("jobstatus");
  if (node && node.dataset.terminal === "0") {
    var key = window.location.pathname.split("/")[2];
    var tick = function () {
      fetch("/run/" + encodeURIComponent(key) + "/fragment/job", { cache: "no-store" })
        .then(function (r) { return r.text(); })
        .then(function (html) {
          var cur = document.getElementById("jobstatus");
          if (!cur) { return; }
          cur.outerHTML = html;
          var next = document.getElementById("jobstatus");
          if (next && next.dataset.terminal === "0") { setTimeout(tick, 5000); }
        })
        .catch(function () { /* a failed poll is left visible as a stale timestamp */ });
    };
    setTimeout(tick, 5000);
  }
})();
