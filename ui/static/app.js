const state = structuredClone(window.DEFAULT_CONFIG);
let activeRunId = null;
let pollTimer = null;
let browserTarget = null;
let browserSelect = "file";
let browserPath = null;
let browserInsert = "";

const $ = (selector, root = document) => root.querySelector(selector);
const $$ = (selector, root = document) => Array.from(root.querySelectorAll(selector));
const NOT_USED_TEXT = "not used";
const NOT_USED_VALUE = "__not_used__";
const LLM_PROVIDERS = {
  openai: {
    baseUrl: "https://api.openai.com/v1/chat/completions",
    model: "gpt-4.1-mini",
    requiresKey: true,
  },
  deepseek: {
    baseUrl: "https://api.deepseek.com/chat/completions",
    model: "deepseek-chat",
    requiresKey: true,
  },
  qwen: {
    baseUrl: "https://dashscope.aliyuncs.com/compatible-mode/v1/chat/completions",
    model: "qwen-plus",
    requiresKey: true,
  },
  moonshot: {
    baseUrl: "https://api.moonshot.cn/v1/chat/completions",
    model: "moonshot-v1-8k",
    requiresKey: true,
  },
  zhipu: {
    baseUrl: "https://open.bigmodel.cn/api/paas/v4/chat/completions",
    model: "glm-4-flash",
    requiresKey: true,
  },
  siliconflow: {
    baseUrl: "https://api.siliconflow.cn/v1/chat/completions",
    model: "Qwen/Qwen2.5-72B-Instruct",
    requiresKey: true,
  },
  openrouter: {
    baseUrl: "https://openrouter.ai/api/v1/chat/completions",
    model: "openai/gpt-4o-mini",
    requiresKey: true,
  },
  ollama: {
    baseUrl: "http://127.0.0.1:11434/v1/chat/completions",
    model: "llama3.1",
    requiresKey: false,
  },
  custom: {
    baseUrl: "",
    model: "",
    requiresKey: false,
  },
};

const HELP_TEXT = {
  programName: "User-facing name for this external program block. It is used in logs to identify the program.",
  executeCommand: "Shell command run inside Command path for each scan point, for example ./TestFunction.py.",
  commandPath: "Working directory where Execute command is run. Relative paths are resolved from the EasyScan root.",
  commandExecutor: "Execution backend. os.system is the default; subprocess.popen is used when a time limit is set and can capture program output.",
  timeLimit: "Maximum runtime in minutes for this program. When set, EasyScan kills the program if it exceeds the limit.",
  cleanOutput: "Whether to delete output files before each program execution.",
  bound: "Additional physical or tabulated bounds. Forms include direct conditions such as x, 0.5, 2.5 or table bounds such as y, x, MAX, bound.txt.",
  fileId: "Numeric ID used to connect variable definitions to an input or output file.",
  filePath: "Path to an input or output file. Relative paths are resolved from the EasyScan root.",
  variableName: "Variable ID used by EasyScan for scanned values, program inputs/outputs, constraints, and plots.",
  variableFileId: "File ID where this variable is read from or written to.",
  variableMethod: "How EasyScan reads or writes the variable: Position, Label, SLHA, Replace, or File.",
  variableArguments: "Method-specific arguments. For Position, use line number and column number, such as 1, 2.",
};

const NUMBER_OF_POINTS_LABELS = {
  RANDOM: {
    label: "Number of points (total number)",
    help: "Total number of random scan points.",
  },
  MCMC: {
    label: "Number of points (total surviving number)",
    help: "Target total number of surviving accepted points in the MCMC scan.",
  },
  BESTFIT: {
    label: "Number of points (max iterations)",
    help: "Maximum differential evolution iterations for the BESTFIT optimizer.",
  },
  EMCEE: {
    label: "Number of points (total samples)",
    help: "Target total saved samples for emcee ensemble MCMC.",
  },
  DYNESTY: {
    label: "Number of points (live points)",
    help: "Number of live points used by dynesty nested sampling.",
  },
  GRID: {
    label: "Number of points (not used)",
    help: "Not used by GRID; use Bins in Input Parameters to set grid density.",
  },
  ONEPOINT: {
    label: "Number of points (not used)",
    help: "Not used by ONEPOINT.",
  },
  ONEPOINTBATCH: {
    label: "Number of points (not used)",
    help: "Not used by ONEPOINTBATCH; points are read from the batch file.",
  },
};

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
      updateConditionalControls();
    });
    el.addEventListener("change", () => {
      setPath(path, el.value);
      updateConditionalControls();
    });
  });
}

function syncRootFields() {
  $$("[data-bind]").forEach((el) => {
    el.value = getPath(el.dataset.bind) ?? "";
  });
}

function replaceState(nextState) {
  Object.keys(state).forEach((key) => delete state[key]);
  Object.assign(state, structuredClone(nextState));
  syncRootFields();
  rerender();
  updateConditionalControls();
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

function helpTip(key) {
  const text = HELP_TEXT[key] || "";
  if (!text) return "";
  return `<span class="help-tip" data-tooltip="${escapeHtml(text)}" tabindex="0" aria-label="Help: ${escapeHtml(key)}">?</span>`;
}

function labelText(text, key) {
  return `<span class="field-label">${escapeHtml(text)} ${helpTip(key)}</span>`;
}

function thText(text, key) {
  return `${escapeHtml(text)} ${helpTip(key)}`;
}

function pathForElement(el) {
  return el?.dataset.field || el?.dataset.bind || "";
}

function restoreControlValue(el) {
  if (!el) return;
  if (el.dataset.originalType) {
    el.type = el.dataset.originalType;
    delete el.dataset.originalType;
  }
  if (el.tagName === "SELECT") {
    $$('option[data-not-used-option="true"]', el).forEach((option) => option.remove());
  }
  const path = pathForElement(el);
  if (path) el.value = getPath(path) ?? "";
}

function showNotUsedValue(el) {
  if (!el) return;
  if (el.tagName === "SELECT") {
    let option = $('option[data-not-used-option="true"]', el);
    if (!option) {
      option = document.createElement("option");
      option.dataset.notUsedOption = "true";
      option.value = NOT_USED_VALUE;
      option.textContent = NOT_USED_TEXT;
      el.prepend(option);
    }
    el.value = NOT_USED_VALUE;
    return;
  }
  if (el.tagName === "INPUT") {
    if (el.type !== "text") {
      el.dataset.originalType ||= el.type;
      el.type = "text";
    }
    el.value = NOT_USED_TEXT;
  } else if (el.tagName === "TEXTAREA") {
    el.value = NOT_USED_TEXT;
  }
}

function setElementDisabled(el, disabled, showNotUsed = false) {
  if (!el) return;
  if (disabled && showNotUsed) {
    showNotUsedValue(el);
  } else {
    restoreControlValue(el);
  }
  el.disabled = disabled;
  el.classList.toggle("not-used-control", disabled && showNotUsed);
}

function setLabelDisabled(label, disabled) {
  if (label) label.classList.toggle("disabled-field", disabled);
}

function setText(selector, text) {
  const el = $(selector);
  if (el) el.textContent = text;
}

function scanMethod() {
  return String(state.scan_method || "").toUpperCase();
}

function updateTopLevelControls() {
  const method = scanMethod();
  const isOnePointBatch = method === "ONEPOINTBATCH";
  const supportsRandomSeed = ["RANDOM", "MCMC", "EMCEE"].includes(method);
  const supportsParallel = ["RANDOM", "GRID", "MCMC", "ONEPOINTBATCH"].includes(method);
  const parallelThreads = Number(state.parallel_threads || 1);
  const usesParallelFolder = supportsParallel && parallelThreads > 1;
  const batchInput = $('[data-bind="batch_file"]');
  const batchButton = $('[data-target="batch_file"]');
  const batchLabel = $("#batchFileLabel");
  const batchLabelText = $("#batchFileLabelText");
  const pointsLabelContainer = $("#numberOfPointsLabel");
  const pointsLabel = $("#numberOfPointsLabelText");
  const pointsHelp = $("#numberOfPointsHelp");
  const pointsInput = $('[data-bind="number_of_points"]');
  const randomSeedInput = $('[data-bind="random_seed"]');
  const parallelThreadsInput = $('[data-bind="parallel_threads"]');
  const parallelFolderInput = $('[data-bind="parallel_folder"]');
  const parallelFolderButton = $('[data-target="parallel_folder"]');
  const pointsText = NUMBER_OF_POINTS_LABELS[method] || NUMBER_OF_POINTS_LABELS.RANDOM;
  const usesNumberOfPoints = ["RANDOM", "BESTFIT", "MCMC", "EMCEE", "DYNESTY"].includes(method);

  setElementDisabled(batchInput, !isOnePointBatch, true);
  setElementDisabled(batchButton, !isOnePointBatch);
  setLabelDisabled(batchLabel, !isOnePointBatch);
  if (batchLabelText) batchLabelText.textContent = isOnePointBatch ? "OnePointBatch file" : "OnePointBatch file (not used)";

  setElementDisabled(pointsInput, !usesNumberOfPoints, true);
  setLabelDisabled(pointsLabelContainer, !usesNumberOfPoints);
  if (pointsLabel) pointsLabel.textContent = pointsText.label;
  if (pointsHelp) pointsHelp.dataset.tooltip = pointsText.help;

  setElementDisabled(randomSeedInput, !supportsRandomSeed, true);
  setLabelDisabled($("#randomSeedLabel"), !supportsRandomSeed);
  setText("#randomSeedLabelText", supportsRandomSeed ? "Random seed" : "Random seed (not used)");

  setElementDisabled(parallelThreadsInput, !supportsParallel, true);
  setLabelDisabled($("#parallelThreadsLabel"), !supportsParallel);
  setText("#parallelThreadsLabelText", supportsParallel ? "Parallel threads" : "Parallel threads (not used)");

  setElementDisabled(parallelFolderInput, !usesParallelFolder, true);
  setElementDisabled(parallelFolderButton, !usesParallelFolder);
  setLabelDisabled($("#parallelFolderLabel"), !usesParallelFolder);
  setText("#parallelFolderLabelText", usesParallelFolder ? "Parallel folder" : "Parallel folder (not used)");
}

function updateInputParameterHeaders() {
  const method = scanMethod();
  const usedByMethod = {
    name: true,
    prior: !["ONEPOINT", "ONEPOINTBATCH"].includes(method),
    value: true,
    minimum: ["RANDOM", "GRID", "BESTFIT", "MCMC", "EMCEE", "DYNESTY"].includes(method),
    maximum: ["RANDOM", "GRID", "BESTFIT", "MCMC", "EMCEE", "DYNESTY"].includes(method),
    bins: method === "GRID",
    interval: ["MCMC", "EMCEE"].includes(method),
    initial: ["MCMC", "EMCEE"].includes(method),
  };
  const labels = {
    name: "Name",
    prior: "Prior",
    value: "Value",
    minimum: "Min",
    maximum: "Max",
    bins: "Bins",
    interval: "Interval",
    initial: "Initial",
  };
  const targets = {
    name: "#inputHeaderName",
    prior: "#inputHeaderPrior",
    value: "#inputHeaderValue",
    minimum: "#inputHeaderMin",
    maximum: "#inputHeaderMax",
    bins: "#inputHeaderBins",
    interval: "#inputHeaderInterval",
    initial: "#inputHeaderInitial",
  };
  Object.entries(targets).forEach(([key, selector]) => {
    setText(selector, `${labels[key]}${usedByMethod[key] ? "" : " (not used)"}`);
    const th = $(selector)?.closest("th");
    setLabelDisabled(th, !usedByMethod[key]);
  });
}

function setFieldDisabled(path, disabled) {
  const el = $(`[data-field="${path}"]`);
  setElementDisabled(el, disabled, true);
  setLabelDisabled(el?.closest("td"), disabled);
}

function updateInputParameterRows() {
  const method = scanMethod();
  state.input_parameters.forEach((parameter, index) => {
    const fixedPrior = String(parameter.prior || "").toUpperCase() === "FIXED";
    const onePointMode = ["ONEPOINT", "ONEPOINTBATCH"].includes(method);
    const enabled = {
      name: true,
      prior: !onePointMode,
      value: onePointMode || fixedPrior,
      minimum: !fixedPrior && ["RANDOM", "GRID", "BESTFIT", "MCMC", "EMCEE", "DYNESTY"].includes(method),
      maximum: !fixedPrior && ["RANDOM", "GRID", "BESTFIT", "MCMC", "EMCEE", "DYNESTY"].includes(method),
      bins: !fixedPrior && method === "GRID",
      interval: !fixedPrior && ["MCMC", "EMCEE"].includes(method),
      initial: !fixedPrior && ["MCMC", "EMCEE"].includes(method),
    };
    setFieldDisabled(`input_parameters.${index}.name`, !enabled.name);
    setFieldDisabled(`input_parameters.${index}.prior`, !enabled.prior);
    setFieldDisabled(`input_parameters.${index}.value`, !enabled.value);
    setFieldDisabled(`input_parameters.${index}.minimum`, !enabled.minimum);
    setFieldDisabled(`input_parameters.${index}.maximum`, !enabled.maximum);
    setFieldDisabled(`input_parameters.${index}.bins`, !enabled.bins);
    setFieldDisabled(`input_parameters.${index}.interval`, !enabled.interval);
    setFieldDisabled(`input_parameters.${index}.initial`, !enabled.initial);
  });
}

function updatePlotRows() {
  state.plots.forEach((plot, index) => {
    const kind = String(plot.kind || "").toUpperCase();
    const enabled = {
      kind: true,
      x: true,
      y: ["SCATTER", "COLOR", "CONTOUR"].includes(kind),
      value: ["COLOR", "CONTOUR"].includes(kind),
      figure_name: true,
    };
    setFieldDisabled(`plots.${index}.kind`, !enabled.kind);
    setFieldDisabled(`plots.${index}.x`, !enabled.x);
    setFieldDisabled(`plots.${index}.y`, !enabled.y);
    setFieldDisabled(`plots.${index}.value`, !enabled.value);
    setFieldDisabled(`plots.${index}.figure_name`, !enabled.figure_name);
  });
}

function updateProgramControls() {
  state.programs.forEach((program, index) => {
    const timeLimitSet = String(program.time_limit || "").trim() !== "";
    const executor = $(`[data-field="programs.${index}.command_executor"]`);
    setElementDisabled(executor, timeLimitSet, true);
    setLabelDisabled(executor?.closest("label"), timeLimitSet);
    setText(`#program${index}ExecutorLabelText`, timeLimitSet ? "Command executor (set by time limit)" : "Command executor");
  });
}

function updateConditionalControls() {
  updateTopLevelControls();
  updateInputParameterHeaders();
  updateInputParameterRows();
  updatePlotRows();
  updateProgramControls();
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
            <tr><th>${thText("ID", "fileId")}</th><th>${thText("Path", "filePath")}</th><th></th></tr>
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
            <tr>
              <th>${thText("Name", "variableName")}</th>
              <th>${thText("File ID", "variableFileId")}</th>
              <th>${thText("Method", "variableMethod")}</th>
              <th>${thText("Arguments", "variableArguments")}</th>
              <th></th>
            </tr>
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
            <label>${labelText("Program name", "programName")}
              ${input("Program name", `programs.${i}.program_name`)}
            </label>
            <label>${labelText("Execute command", "executeCommand")}
              ${input("Execute command", `programs.${i}.execute_command`)}
            </label>
            <label>${labelText("Command path", "commandPath")}
              <div class="path-row">
                ${input("Command path", `programs.${i}.command_path`)}
                <button class="icon-button browse" data-target="programs.${i}.command_path" data-select="dir" title="Choose folder">...</button>
              </div>
            </label>
            <label><span class="field-label"><span id="program${i}ExecutorLabelText">Command executor</span> ${helpTip("commandExecutor")}</span>
              ${select(`programs.${i}.command_executor`, ["", "os.system", "subprocess.popen"])}
            </label>
            <label>${labelText("Time limit in minutes", "timeLimit")}
              ${input("Time limit", `programs.${i}.time_limit`)}
            </label>
            <label>${labelText("Clean output file", "cleanOutput")}
              ${select(`programs.${i}.clean_output`, ["", "Yes", "No"])}
            </label>
          </div>
          ${renderFileRows(i, "input_files", "Input Files", "file")}
          ${renderVariableRows(i, "input_variables", "Input Variables")}
          ${renderFileRows(i, "output_files", "Output Files", "file")}
          ${renderVariableRows(i, "output_variables", "Output Variables")}
          <div class="subsection">
            <label>${labelText("Bound", "bound")}
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
    el.addEventListener("input", () => {
      setPath(el.dataset.field, el.value);
      updateConditionalControls();
    });
    el.addEventListener("change", () => {
      setPath(el.dataset.field, el.value);
      updateConditionalControls();
    });
  });
}

function rerender() {
  renderAll();
  refreshBindings();
  updateConditionalControls();
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

function setConfigFileStatus(message, kind = "") {
  const el = $("#configFileStatus");
  if (!el) return;
  el.textContent = message;
  el.className = `config-status ${kind}`.trim();
}

function setAiConfigStatus(message, kind = "") {
  const el = $("#aiConfigStatus");
  if (!el) return;
  el.textContent = message;
  el.className = `config-status ${kind}`.trim();
}

function allProviderBaseUrls() {
  return Object.values(LLM_PROVIDERS)
    .map((item) => item.baseUrl)
    .filter(Boolean);
}

function applyLlmProviderDefaults(force = false) {
  const providerEl = $("#aiProvider");
  const modelEl = $("#aiModel");
  const apiKeyEl = $("#aiApiKey");
  const baseUrlEl = $("#aiBaseUrl");
  if (!providerEl || !modelEl || !apiKeyEl || !baseUrlEl) return;
  const defaults = LLM_PROVIDERS[providerEl.value] || LLM_PROVIDERS.custom;
  if (force || !modelEl.value.trim()) {
    modelEl.value = defaults.model;
  }
  if (force || !baseUrlEl.value.trim() || allProviderBaseUrls().includes(baseUrlEl.value.trim())) {
    baseUrlEl.value = defaults.baseUrl;
  }
  apiKeyEl.disabled = !defaults.requiresKey;
  apiKeyEl.placeholder = defaults.requiresKey ? "Required" : "Not required";
  if (!defaults.requiresKey) {
    apiKeyEl.value = "";
  }
}

function detailMessage(detail, fallback) {
  if (!detail) return fallback;
  if (typeof detail === "string") return detail;
  if (detail.message) return detail.message;
  if (detail.check_text) return detail.check_text;
  try {
    return JSON.stringify(detail);
  } catch {
    return fallback;
  }
}

function showCheckReport(text, ok) {
  const box = $("#checkBox");
  if (!box || !text) return;
  box.classList.remove("hidden");
  box.textContent = text;
  box.classList.toggle("failed", !ok);
  box.classList.toggle("completed", Boolean(ok));
}

async function importConfig() {
  setConfigFileStatus("Loading...");
  const response = await fetch("/api/config/import/open", { method: "POST" });
  const payload = await response.json();
  if (!response.ok) {
    setConfigFileStatus(payload.detail || "Failed to load INI.", response.status === 400 ? "" : "failed");
    return;
  }
  replaceState(payload);
  setConfigFileStatus(`Loaded ${payload.config_file_path || "INI file"}.`, "completed");
}

async function exportConfig() {
  setConfigFileStatus("Exporting...");
  const response = await fetch("/api/config/export/save", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(state),
  });
  const payload = await response.json().catch(() => ({}));
  if (!response.ok) {
    setConfigFileStatus(payload.detail || "Failed to export INI.", response.status === 400 ? "" : "failed");
    return;
  }
  if (payload.path) state.config_file_path = payload.path;
  setConfigFileStatus(`Saved ${payload.path || "INI file"}.`, "completed");
}

async function checkConfig() {
  setConfigFileStatus("Checking...");
  const response = await fetch("/api/config/check", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(state),
  });
  const payload = await response.json().catch(() => ({}));
  const box = $("#checkBox");
  if (box) {
    box.classList.remove("hidden");
    box.textContent = payload.text || payload.detail || "Check failed.";
    box.classList.toggle("failed", !payload.ok);
    box.classList.toggle("completed", Boolean(payload.ok));
  }
  if (!response.ok) {
    setConfigFileStatus(payload.detail || "Check failed.", "failed");
    return;
  }
  setConfigFileStatus(payload.ok ? "Config check passed." : "Config check found issues.", payload.ok ? "completed" : "failed");
}

async function generateConfigWithAI() {
  const button = $("#aiGenerateConfigBtn");
  const prompt = $("#aiPrompt")?.value.trim() || "";
  if (!prompt) {
    setAiConfigStatus("Enter a natural-language request first.", "failed");
    return;
  }

  if (button) button.disabled = true;
  setAiConfigStatus("Generating...");
  setConfigFileStatus("Generating INI with AI...");
  const response = await fetch("/api/config/ai/generate", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      provider: $("#aiProvider")?.value || "openai",
      model: $("#aiModel")?.value.trim() || "",
      api_key: $("#aiApiKey")?.value.trim() || "",
      base_url: $("#aiBaseUrl")?.value.trim() || "",
      prompt,
      current_config: state,
    }),
  });
  const payload = await response.json().catch(() => ({}));
  if (!response.ok) {
    const detail = payload.detail || {};
    if (detail.check_text) {
      showCheckReport(detail.check_text, false);
    } else if (detail.ini) {
      showCheckReport(detail.ini, false);
    }
    const message = detailMessage(detail, "Failed to generate INI.");
    setAiConfigStatus(message, "failed");
    setConfigFileStatus(message, "failed");
    if (button) button.disabled = false;
    return;
  }

  if (payload.config) {
    const existingPath = state.config_file_path;
    if (existingPath && !payload.config.config_file_path) {
      payload.config.config_file_path = existingPath;
    }
    replaceState(payload.config);
  }
  showCheckReport(payload.check_text || "Config check passed.", payload.check?.ok ?? true);
  setAiConfigStatus("Generated and loaded into the UI.", "completed");
  setConfigFileStatus("AI generated config loaded into the UI.", "completed");
  if (button) button.disabled = false;
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

async function openBrowser(target, selectKind, insertMode = "", startPath = null) {
  browserTarget = target;
  browserSelect = selectKind;
  browserInsert = insertMode;
  browserPath = startPath || getPath(target) || null;
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

async function choosePath(path) {
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
    await choosePath(button.dataset.pickPath);
  }
});

$("#addInputParam").addEventListener("click", () => addRow("input_parameters"));
$("#addProgram").addEventListener("click", () => addRow("programs"));
$("#addGaussian").addEventListener("click", () => addRow("gaussian_constraints"));
$("#addPlot").addEventListener("click", () => addRow("plots"));
$("#importConfigBtn").addEventListener("click", importConfig);
$("#exportConfigBtn").addEventListener("click", exportConfig);
$("#checkConfigBtn").addEventListener("click", checkConfig);
$("#aiProvider").addEventListener("change", () => applyLlmProviderDefaults(true));
$("#aiGenerateConfigBtn").addEventListener("click", generateConfigWithAI);
$("#runBtn").addEventListener("click", startRun);
$("#stopBtn").addEventListener("click", stopRun);
$("#closeBrowser").addEventListener("click", () => $("#browserModal").classList.add("hidden"));
$("#upBtn").addEventListener("click", () => loadBrowser($("#upBtn").dataset.path));
$("#rootSelect").addEventListener("change", (event) => loadBrowser(event.target.value));
$("#selectCurrentBtn").addEventListener("click", () => choosePath(browserPath));

bindRootFields();
applyLlmProviderDefaults(true);
rerender();
updateConditionalControls();
