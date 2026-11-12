'use strict';
const assert = require('node:assert/strict');
const {runtime, dictionaries} = require('./u4_web_contract.cjs');

async function labels() {
  const {api, document} = await runtime();
  const expected = new Map([
    ['query-csv', 's.063'], ['query-smiles-col', 's.052'], ['query-id-col', 's.065'],
    ['new-activity-col', 'u8.activityColumn'], ['new-unit', 's.054'], ['new-unit-col', 's.055'],
    ['new-relation-col', 's.058'], ['new-document-col', 's.056'], ['new-label-col', 'u8.labelColumn'],
    ['new-class-filter', 'u8.categoryFilter'], ['new-label-filter', 'u8.labelFilter'],
    ...['smiles','id','label','mode','value','unit','class','batch','note'].map(name => ['single-'+name,'u9b.'+name]),
  ]);
  assert.equal(document.querySelectorAll('[data-help]').length, 80);
  assert.equal(document.querySelectorAll('dialog [data-help]').length, expected.size + 1);
  for (const language of ['en', 'zh', 'en']) {
    api.setLanguage(language);
    for (const [control, key] of expected) {
      const label = document.querySelector(`label[for="${control}"]`);
      assert.ok(label, control);
      assert.equal(label.childNodes[0].nodeType, 3, `${control}: visible name precedes help`);
      assert.equal(label.childNodes[0].textContent, dictionaries[language][key]);
      const button = label.querySelector('[data-help]');
      assert.equal(label.children[0], button, `${control}: question mark follows name`);
      assert.equal(button.textContent, '?');
      assert.equal(button.getAttribute('aria-label'), api.t('u.help', {name: dictionaries[language][key]}));
      assert.ok(label.contains(document.getElementById(control)), `${control}: field preserved`);
    }
    assert.equal(document.querySelector('label[for="query-csv"]').childNodes[0].textContent,
      language === 'en' ? 'CSV or TSV file' : 'CSV 或 TSV 文件');
  }
}

async function dialogHelp() {
  const {api, document, window} = await runtime();
  const tooltip = document.getElementById('parameter-tooltip');
  let checked = 0;
  for (const dialog of document.querySelectorAll('dialog')) {
    const buttons = dialog.querySelectorAll('[data-help]');
    if (!buttons.length) continue;
    dialog.showModal();
    for (const button of buttons) {
      checked++;
      document.activeElement = null;
      button.dispatchEvent({type: 'pointerenter', pointerType: 'mouse'});
      assert.equal(tooltip.hidden, false);
      assert.equal(tooltip.parentNode, dialog, `${button.dataset.helpControl}: shares modal top layer`);
      assert.equal(document.querySelectorAll('[role=tooltip]').length, 1, 'Reparent the shared tooltip');
      assert.equal(tooltip.textContent, api.t(button.dataset.help, {activity: 'pChEMBL'}));
      button.dispatchEvent({type: 'pointerleave', relatedTarget: tooltip});
      assert.equal(tooltip.hidden, false, 'Moving into tooltip keeps it open');
      tooltip.dispatchEvent({type: 'pointerleave', relatedTarget: dialog});
      assert.equal(tooltip.hidden, true);
      button.focus();
      assert.equal(tooltip.hidden, false, 'Keyboard focus opens help');
      const firstEscape = {type: 'keydown', key: 'Escape'};
      document.dispatchEvent(firstEscape);
      assert.equal(tooltip.hidden, true, 'Escape closes help');
      assert.equal(firstEscape.defaultPrevented, true, 'First Escape preserves the open dialog');
      const nextEscape = {type: 'keydown', key: 'Escape'};
      document.dispatchEvent(nextEscape);
      assert.equal(!!nextEscape.defaultPrevented, false, 'Next Escape can close the dialog normally');
      button.focus();
      button.dispatchEvent({type: 'blur', relatedTarget: dialog});
      assert.equal(tooltip.hidden, true, 'Leaving keyboard focus closes help');
      button.click();
      assert.equal(tooltip.hidden, false);
      button.click();
      assert.equal(tooltip.hidden, true, 'Pinned help toggles closed');
      for (const width of [360, 1280]) {
        window.innerWidth = width;
        button.rect = {left: width - 20, right: width - 4, top: 780, bottom: 796, width: 16, height: 16};
        api.showHelp(button);
        assert.ok(parseFloat(tooltip.style.left) >= 0);
        assert.ok(parseFloat(tooltip.style.left) + tooltip.offsetWidth <= width);
        assert.ok(parseFloat(tooltip.style.top) + tooltip.offsetHeight <= window.innerHeight);
      }
      dialog.close();
      assert.equal(tooltip.hidden, true, 'Dialog closing clears its help');
      dialog.showModal();
    }
    dialog.close();
  }
  assert.equal(checked, 21, 'Every question mark in every dialog checked');
  const outside = document.querySelector('[data-help-control="input-format"]');
  api.showHelp(outside);
  assert.equal(tooltip.parentNode, document.body, 'Help returns to body outside dialogs');
  assert.equal(document.querySelectorAll('[role=tooltip]').length, 1);
  const outsideEscape = {type: 'keydown', key: 'Escape'};
  document.dispatchEvent(outsideEscape);
  assert.equal(tooltip.hidden, true);
  assert.equal(!!outsideEscape.defaultPrevented, false, 'Non-dialog help keeps the U4 Escape behavior');
}

const scenarios = {labels, dialogHelp};
(async () => {
  const scenario = process.argv[3];
  assert.ok(scenario in scenarios, scenario);
  await scenarios[scenario]();
  console.log(`U8b ${scenario}: passed`);
})().catch(error => { console.error(error); process.exitCode = 1; });
