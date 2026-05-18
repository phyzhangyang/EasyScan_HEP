const state = structuredClone(window.DEFAULT_CONFIG);
let activeRunId = null;
let pollTimer = null;
let browserTarget = null;
let browserSelect = "file";
let browserPath = null;
let browserInsert = "";

const $ = (selector, root = document) => root.querySelector(selector);
const $$ = (selector, root = document) => Array.from(root.querySelectorAll(selector));

function getPath(path) {
  return path.split(".").reduce((obj, part) => obj?.[part], state);
}

function setPath(path, value) {
  const parts = path.split(".");
  let obj = state;
  for (let i = 0; i < parts.length - 1; i += 1) {
    obj = obj[parts[i]];
  }
  obj[parts.at(-1)] = value;
}

function bindRootFields() {
  $$("[data-bind]").forEach((el) => {
    const path = el.dataset.bind;
    el.value = getPath(path) ?? "";
    el.addEventListener("input", () => {
      setPath(path, el.value);
    });
    el.addEventListener("change", () => {
      setPath(path, el.value);
    });
  });
}

function input(name, path, type = "text", placeholder = "") {
  const value = getPath(path) ?? "";
  return `<input aria-label="${name}" data-field="${path}" type="${type}" value="${escapeHtml(value)}" placeholder="${placeholder}">`;
}

function select(path, options) {
  const value = getPath(path) ?? "";
  const opts = options
    .map((item) => `<option ${item === value ? "selected" : ""}>${item}</option>`)
    .join("");
  return `<select data-field="${path}">${opts}</select>`;
}

function escapeHtml(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;");
}

function renderInputParameters() {
  $("#inputParamsBody").innerHTML = state.input_parameters
    .map(
      (_, i) => `
        <tr>
          <td>${input("Name", `input_parameters.${i}.name`)}</td>
          <td>${select(`input_parameters.${i}.prior`, ["flat", "log", "Fixed"])}</td>
          <td>${input("Value", `input_parameters.${i}.value`)}</td>
          <td>${input("Min", `input_parameters.${i}.minimum`)}</td>
          <td>${input("Max", `input_parameters.${i}.maximum`)}</td>
          <td>${input("Bins", `input_parameters.${i}.bins`)}</td>
          <td>${input("Interval", `input_parameters.${i}.interval`)}</td>
          <td>${input("Initial", `input_parameters.${i}.initial`)}</td>
          <td><button class="danger small" data-remove="input_parameters.${i}">Remove</button></td>
        </tr>
      `
    )
    .join("");
}

function renderFileRows(programIndex, fieldName, label, selectKind) {
  const rows = state.programs[programIndex][fieldName];
  return `
    <div class="subsection">
      <div class="subsection-title">
        <h3>${label}</h3>
        <button class="secondary small" data-add="${fieldName}" data-program="${programIndex}">Add</button>
      </div>
      <div class="table-wrap">
        <table>
          <thead>
            <tr><th>ID</th><th>Path</th><th></th></tr>
          </thead>
          <tbody>
            ${rows
              .map(
                (_, rowIndex) => `
                  <tr>
                    <td>${input("ID", `programs.${programIndex}.${fieldName}.${rowIndex}.file_id`)}</td>
                    <td>
                      <div class="path-row">
                        ${input("Path", `programs.${programIndex}.${fieldName}.${rowIndex}.path`)}
                        <button class="icon-button browse" data-target="programs.${programIndex}.${fieldName}.${rowIndex}.path" data-select="${selectKind}" title="Choose path">...</button>
                      </div>
                    </td>
                    <td><button class="danger small" data-remove="programs.${programIndex}.${fieldName}.${rowIndex}">Remove</button></td>
                  </tr>
                `
              )
              .join("")}
          </tbody>
        </table>
      </div>
    </div>
  `;
}

function renderVariableRows(programIndex, fieldName, label) {
  const rows = state.programs[programIndex][fieldName];
  return `
    <div class="subsection">
      <div class="subsection-title">
        <h3>${label}</h3>
        <button class="secondary small" data-add="${fieldName}" data-program="${programIndex}">Add</button>
      </div>
      <div class="table-wrap">
        <table>
          <thead>
            <tr><th>Name</th><th>File ID</th><th>Method</th><th>Arguments</th><th></th></tr>
          </thead>
          <tbody>
            ${rows
              .map(
                (_, rowIndex) => `
                  <tr>
                    <td>${input("Name", `programs.${programIndex}.${fieldName}.${rowIndex}.name`)}</td>
                    <td>${input("File ID", `programs.${programIndex}.${fieldName}.${rowIndex}.file_id`)}</td>
                    <td>${select(`programs.${programIndex}.${fieldName}.${rowIndex}.method`, ["Position", "Label", "SLHA", "Replace", "File"])}</td>
                    <td>${input("Arguments", `programs.${programIndex}.${fieldName}.${rowIndex}.arguments`, "text", "Example: 1, 1")}</td>
                    <td><button class="danger small" data-remove="programs.${programIndex}.${fieldName}.${rowIndex}">Remove</button></td>
                  </tr>
                `
              )
              .join("")}
          </tbody>
        </table>
      </div>
    </div>
  `;
}

function renderPrograms() {
  $("#programs").innerHTML = state.programs
    .map(
      (_, i) => `
        <div class="program-card">
          <div class="program-header">
            <h3>Program ${i + 1}</h3>
            <button class="danger small" data-remove="programs.${i}">Remove Program</button>
          </div>
          <div class="grid">
            <label>Program name
              ${input("Program name", `programs.${i}.program_name`)}
            </label>
            <label>Execute command
              ${input("Execute command", `programs.${i}.execute_command`)}
            </label>
            <label>Command path
              <div class="path-row">
                ${input("Command path", `programs.${i}.command_path`)}
                <button class="icon-button browse" data-target="programs.${i}.command_path" data-select="dir" title="Choose folder">...</button>
              </div>
            </label>
            <label>Command executor
              ${select(`programs.${i}.command_executor`, ["", "os.system", "subprocess.popen"])}
            </label>
            <label>Time limit in minutes
              ${input("Time limit", `programs.${i}.time_limit`)}
            </label>
            <label>Clean output file
              ${select(`programs.${i}.clean_output`, ["", "Yes", "No"])}
            </label>
          </div>
          ${renderFileRows(i, "input_files", "Input Files", "file")}
          ${renderVariableRows(i, "input_variables", "Input Variables")}
          ${renderFileRows(i, "output_files", "Output Files", "file")}
          ${renderVariableRows(i, "output_variables", "Output Variables")}
          <div class="subsection">
            <label>Bound
              <div class="path-row">
                <textarea data-field="programs.${i}.bounds" placeholder="Example: x, y, MAX, bounds.dat">${escapeHtml(state.programs[i].bounds)}</textarea>
                <button class="icon-button browse" data-target="programs.${i}.bounds" data-select="file" data-insert="bound-file" title="Choose bound file">...</button>
              </div>
            </label>
          </div>
        </div>
      `
    )
    .join("");
}

function renderGaussian() {
  $("#gaussianBody").innerHTML = state.gaussian_constraints
    .map(
      (_, i) => `
        <tr>
          <td>${input("Variable", `gaussian_constraints.${i}.variable`)}</td>
          <td>${input("Mean", `gaussian_constraints.${i}.mean`)}</td>
          <td>${input("Deviation", `gaussian_constraints.${i}.deviation`)}</td>
          <td>${select(`gaussian_constraints.${i}.kind`, ["symm", "upper", "lower"])}</td>
          <td>${input("Name", `gaussian_constraints.${i}.label`)}</td>
          <td><button class="danger small" data-remove="gaussian_constraints.${i}">Remove</button></td>
        </tr>
      `
    )
    .join("");
}

function renderPlots() {
  $("#plotsBody").innerHTML = state.plots
    .map(
      (_, i) => `
        <tr>
          <td>${select(`plots.${i}.kind`, ["Histogram", "Scatter", "Color", "Contour"])}</td>
          <td>${input("X", `plots.${i}.x`)}</td>
          <td>${input("Y", `plots.${i}.y`)}</td>
          <td>${input("Value", `plots.${i}.value`)}</td>
          <td>${input("Figure", `plots.${i}.figure_name`)}</td>
          <td><button class="danger small" data-remove="plots.${i}">Remove</button></td>
        </tr>
      `
    )
    .join("");
}

function renderAll() {
  renderInputParameters();
  renderPrograms();
  renderGaussian();
  renderPlots();
}

function refreshBindings() {
  $$("[data-field]").forEach((el) => {
    el.addEventListener("input", () => setPath(el.dataset.field, el.value));
    el.addEventListener("change", () => setPath(el.dataset.field, el.value));
  });
}

function rerender() {
  renderAll();
  refreshBindings();
}

function removePath(path) {
  const parts = path.split(".");
  const index = Number(parts.pop());
  const arr = getPath(parts.join("."));
  if (Array.isArray(arr) && arr.length > 1) {
    arr.splice(index, 1);
  }
  rerender();
}

function addRow(kind, programIndex = null) {
  if (kind === "input_parameters") {
    state.input_parameters.push({ name: "new_parameter", prior: "flat", value: "", minimum: "0", maximum: "1", bins: "10", interval: "10", initial: "0.5" });
  } else if (kind === "programs") {
    state.programs.push(structuredClone(state.programs[0]));
  } else if (kind === "gaussian_constraints") {
    state.gaussian_constraints.push({ variable: "", mean: "", deviation: "", kind: "symm", label: "" });
  } else if (kind === "plots") {
    state.plots.push({ kind: "Color", x: "", y: "", value: "", figure_name: "" });
  } else if (programIndex !== null) {
    const program = state.programs[programIndex];
    if (kind === "input_files" || kind === "output_files") {
      program[kind].push({ file_id: String(program[kind].length + 1), path: "" });
    } else {
      program[kind].push({ name: "", file_id: "1", method: "Position", arguments: "" });
    }
  }
  rerender();
}

async function previewConfig() {
  const response = await fetch("/api/config", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(state),
  });
  const payload = await response.json();
  $("#configPreview").textContent = payload.ini || JSON.stringify(payload, null, 2);
}

async function startRun() {
  $("#runBtn").disabled = true;
  $("#stopBtn").disabled = false;
  $("#logBox").textContent = "";
  $("#results").innerHTML = "";
  setStatus("queued");
  const response = await fetch("/api/runs", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(state),
  });
  const payload = await response.json();
  if (!response.ok) {
    $("#logBox").textContent = payload.detail || "Failed to start run.";
    setStatus("failed");
    $("#runBtn").disabled = false;
    $("#stopBtn").disabled = true;
    return;
  }
  activeRunId = payload.run_id;
  pollRun();
  pollTimer = setInterval(pollRun, 1200);
}

async function pollRun() {
  if (!activeRunId) return;
  const response = await fetch(`/api/runs/${activeRunId}`);
  const payload = await response.json();
  if (!response.ok) return;
  setStatus(payload.status);
  $("#logBox").textContent = payload.logs || "";
  $("#logBox").scrollTop = $("#logBox").scrollHeight;
  renderResults(payload.results);
  if (["completed", "failed", "stopped"].includes(payload.status)) {
    clearInterval(pollTimer);
    $("#runBtn").disabled = false;
    $("#stopBtn").disabled = true;
  }
}

async function stopRun() {
  if (!activeRunId) return;
  await fetch(`/api/runs/${activeRunId}/stop`, { method: "POST" });
  setStatus("stopped");
}

function setStatus(status) {
  const el = $("#runStatus");
  el.textContent = status.charAt(0).toUpperCase() + status.slice(1);
  el.className = `status ${status}`;
}

function renderResults(results) {
  if (!results) return;
  const files = results.files || [];
  const plots = results.plots || [];
  $("#results").innerHTML = `
    <div><strong>Result directory:</strong> ${escapeHtml(results.result_dir || "")}</div>
    <div class="result-list">
      ${files.map((item) => `<a href="${item.url}" target="_blank">${escapeHtml(item.name)}</a>`).join("")}
    </div>
    <div class="plot-grid">
      ${plots
        .map(
          (item) => `
            <a class="plot-card" href="${item.url}" target="_blank">
              <img src="${item.url}" alt="${escapeHtml(item.name)}">
              <div>${escapeHtml(item.name)}</div>
            </a>
          `
        )
        .join("")}
    </div>
  `;
}

async function openBrowser(target, selectKind, insertMode = "") {
  browserTarget = target;
  browserSelect = selectKind;
  browserInsert = insertMode;
  browserPath = getPath(target) || null;
  $("#browserModal").classList.remove("hidden");
  await loadBrowser(browserPath);
}

async function loadBrowser(path = null) {
  const params = new URLSearchParams({ select: browserSelect });
  if (path) params.set("path", path);
  const response = await fetch(`/api/browse?${params.toString()}`);
  const payload = await response.json();
  browserPath = payload.current;
  $("#currentPath").textContent = payload.current;
  $("#rootSelect").innerHTML = payload.roots.map((root) => `<option value="${escapeHtml(root.path)}">${escapeHtml(root.name)}</option>`).join("");
  $("#rootSelect").value = payload.roots.find((root) => payload.current.startsWith(root.path))?.path || payload.roots[0]?.path || "";
  $("#upBtn").dataset.path = payload.parent;
  $("#selectCurrentBtn").disabled = browserSelect !== "dir";
  $("#browserEntries").innerHTML = payload.entries
    .map(
      (entry) => `
        <div class="browser-entry">
          <div class="entry-name">${entry.is_dir ? "[Folder]" : "[File]"} ${escapeHtml(entry.name)}</div>
          <button class="secondary small" data-open-path="${escapeHtml(entry.path)}" ${entry.is_dir ? "" : "disabled"}>Open</button>
          <button class="primary small" data-pick-path="${escapeHtml(entry.path)}" ${entry.selectable ? "" : "disabled"}>Choose</button>
        </div>
      `
    )
    .join("");
}

function choosePath(path) {
  if (!browserTarget) return;
  if (browserInsert === "bound-file") {
    const current = String(getPath(browserTarget) || "").trimEnd();
    const separator = current && !current.endsWith(",") ? ", " : "";
    setPath(browserTarget, `${current}${separator}${path}`);
  } else {
    setPath(browserTarget, path);
  }
  const bound = $(`[data-bind="${browserTarget}"]`) || $(`[data-field="${browserTarget}"]`);
  if (bound) bound.value = getPath(browserTarget);
  $("#browserModal").classList.add("hidden");
  browserTarget = null;
  browserInsert = "";
}

document.addEventListener("click", async (event) => {
  const button = event.target.closest("button");
  if (!button) return;
  if (button.matches(".browse")) {
    await openBrowser(button.dataset.target, button.dataset.select || "file", button.dataset.insert || "");
  } else if (button.dataset.remove) {
    removePath(button.dataset.remove);
  } else if (button.dataset.add) {
    addRow(button.dataset.add, button.dataset.program ? Number(button.dataset.program) : null);
  } else if (button.dataset.openPath) {
    await loadBrowser(button.dataset.openPath);
  } else if (button.dataset.pickPath) {
    choosePath(button.dataset.pickPath);
  }
});

$("#addInputParam").addEventListener("click", () => addRow("input_parameters"));
$("#addProgram").addEventListener("click", () => addRow("programs"));
$("#addGaussian").addEventListener("click", () => addRow("gaussian_constraints"));
$("#addPlot").addEventListener("click", () => addRow("plots"));
$("#previewBtn").addEventListener("click", previewConfig);
$("#runBtn").addEventListener("click", startRun);
$("#stopBtn").addEventListener("click", stopRun);
$("#copyConfig").addEventListener("click", () => navigator.clipboard.writeText($("#configPreview").textContent));
$("#closeBrowser").addEventListener("click", () => $("#browserModal").classList.add("hidden"));
$("#upBtn").addEventListener("click", () => loadBrowser($("#upBtn").dataset.path));
$("#rootSelect").addEventListener("change", (event) => loadBrowser(event.target.value));
$("#selectCurrentBtn").addEventListener("click", () => choosePath(browserPath));

bindRootFields();
rerender();
previewConfig();
