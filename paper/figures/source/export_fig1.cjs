// Export Fig. 1 (v0.7) from its SVG to PDF and a 600-dpi PNG with the pre-installed Chromium.
const {chromium} = require('playwright');
const fs = require('node:fs');
(async () => {
  // Usage: node export_fig1.cjs <folder with Figure_1_workflow.svg> [output folder]
  const dir = process.argv[2], out = process.argv[3] || dir;
  fs.mkdirSync(out, {recursive: true});
  const svg = fs.readFileSync(dir + '/Figure_1_workflow.svg', 'utf8');
  const heightMm = Number(svg.match(/height="([\d.]+)mm"/)[1]);
  const html = `<!doctype html><html><head><meta charset="utf-8"><style>@page{size:170mm ${heightMm}mm;margin:0}html,body{margin:0;padding:0}svg{display:block;width:170mm;height:${heightMm}mm}</style></head><body>${svg.replace(/<\?xml[^>]*>/, '')}</body></html>`;
  const browser = await chromium.launch({headless: true});
  const cssWidth = 170 / 25.4 * 96, cssHeight = heightMm / 25.4 * 96;
  const scale = (170 / 25.4 * 600) / cssWidth;
  const page = await browser.newPage({viewport: {width: Math.ceil(cssWidth), height: Math.ceil(cssHeight)}, deviceScaleFactor: scale});
  await page.setContent(html);
  await page.pdf({path: out + '/Figure_1.pdf', width: '170mm', height: `${heightMm}mm`, printBackground: true});
  await page.locator('svg').screenshot({path: out + '/Figure_1.png', type: 'png'});
  await browser.close();
  console.log('scale', scale, 'height', heightMm);
})().catch(e => { console.error(e); process.exit(1); });
