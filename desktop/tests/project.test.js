'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs/promises');
const path = require('node:path');
const os = require('node:os');
const { zipDirectory, extractZip } = require('../src/core/archive');
const { MANIFEST, writeProjectManifest, readProject, openProject } = require('../src/core/project');

const MESH = ['CM2D 1','VERTICES 4','0 0 0','1 1 0','2 1 1','3 0 1','EDGES 4',
  '0 0 1 0 -1 2','1 1 2 0 -1 2','2 2 3 0 -1 2','3 3 0 0 -1 2','CELLS 1',
  '0 0 0 1 4 0 1 2 3 4 0 1 2 3','AUDIT 0',''].join('\n');
const REQUEST = { case:'channel', nu:.1, speed:1, maxIterations:100, convection:'upwind',
  pressurePreconditioner:'ic0', mode:'transient', dt:.01, steps:2, diffusivity:.1, initial:300,
  source:0, scalarConvection:'upwind', boundaries:Object.fromEntries(
    ['wall','inlet','outlet','top','bottom'].map(key => [key,{kind:'value',value:300,inflowValue:300}])) };

async function fixture(t) {
  const temp = await fs.mkdtemp(path.join(os.tmpdir(), 'project roundtrip '));
  t.after(() => fs.rm(temp, { recursive:true, force:true }));
  const root = path.join(temp, '原始项目'), session = path.join(temp, '新会话');
  await fs.mkdir(root); await fs.mkdir(session);
  const prefix = path.join(root, 'rectangle');
  await fs.writeFile(prefix+'.xy', '0 0\n1 0\n1 1\n0 1\n');
  await fs.writeFile(prefix+'.solver.cm2d', MESH);
  await fs.mkdir(prefix+'-openfoam');
  const current = { outputDirectory:root, prefix, cm2dPath:prefix+'.solver.cm2d',
    job:{geometryPath:'/old/computer/input.dxf',outputDirectory:root,sourceUnits:'mm',method:'cutcell',fluidRegion:'interior'},
    result:{counts:{cells:1},openFoam:{path:prefix+'-openfoam'}},levelBasis:'size',incomplete:null };
  return {temp,root,session,prefix,current};
}

test('project moves to a new session without original geometry, mesh or checkpoint paths', async t => {
  const f = await fixture(t);
  const external = path.join(f.temp,'外部.thermal.checkpoint');
  const checkpoint = ['CARTMESH2D_THERMAL_CHECKPOINT 1','COUPLING new-time-flux-Euler-v1','SCALAR 1 300',
    'FLOW','CARTMESH2D_FLOW_CHECKPOINT 2','CONFIG "channel" 0.1 1 upwind symmetric 0 reject',
    'TIME 0.02','FLUX 4 0 0 0 0',''].join('\r\n');
  await fs.writeFile(external,checkpoint);
  f.current.thermalRestart = {path:external,metadata:{time:.02,request:REQUEST}};
  const ui = {inputs:{flowDt:'0.025',thermalResume:true,thermalInitial:'300'},regions:[]};
  const manifest = await writeProjectManifest(f.current, ui);
  assert.equal(manifest.state.job.sourceUnits,'m');
  assert.equal(manifest.state.thermalRestart.path,'project://project-inputs/thermal.checkpoint');
  assert.equal(f.current.thermalRestart.path,external,'saving does not rebind the live accepted state');
  const archive = path.join(f.temp,'portable.zip');
  await zipDirectory(f.root,archive);
  await fs.rm(f.root,{recursive:true}); await fs.unlink(external);
  const loaded = await openProject(archive,f.session);
  assert.equal(loaded.mesh.cells.length,1);
  assert.equal(loaded.thermalRestart.metadata.time,.02);
  assert.equal(loaded.projectUi.inputs.flowDt,'0.025');
  assert.equal(await fs.readFile(loaded.thermalRestart.path,'utf8'),checkpoint);
  assert.ok(loaded.job.geometryPath.startsWith(f.session+path.sep));
  assert.ok(loaded.result.openFoam.path.startsWith(f.session+path.sep));
  assert.equal(loaded.job.sourceUnits,'m','converted XY cannot be rescaled as the original DXF');
  await writeProjectManifest(loaded,null);
  const second = await readProject(loaded.outputDirectory);
  assert.deepEqual(second.projectUi,ui,'saving with the preview released retains the previous controls');
});

test('corrupt input and missing files reject before publishing a restored project', async t => {
  const f = await fixture(t);
  await writeProjectManifest(f.current,{});
  await fs.appendFile(f.current.cm2dPath,'\n');
  await assert.rejects(readProject(f.root),/文件损坏或被更改/);
  const zip = path.join(f.temp,'bad.zip');
  await zipDirectory(f.root,zip);
  await assert.rejects(openProject(zip,f.session),/文件损坏或被更改/);
  assert.deepEqual(await fs.readdir(f.session),[],'failed extraction is not left as a selectable project');
  assert.equal(await fs.readFile(f.current.cm2dPath,'utf8'),MESH+'\n','original data is untouched');
});

test('old exports and escaping manifest references fail with a useful reason', async t => {
  const f = await fixture(t);
  const old = path.join(f.temp,'old.zip');
  await zipDirectory(f.root,old);
  await assert.rejects(openProject(old,f.session),/没有项目清单/);
  const manifest = await writeProjectManifest(f.current,{});
  manifest.state.job.geometryPath='project://../external.xy';
  await fs.writeFile(path.join(f.root,MANIFEST),JSON.stringify(manifest));
  await assert.rejects(readProject(f.root),/越出项目目录/);
});

test('project open cancellation preserves the source archive', async t => {
  const f = await fixture(t);
  await writeProjectManifest(f.current,{});
  const zip = path.join(f.temp,'cancel.zip');
  await zipDirectory(f.root,zip);
  const before=await fs.readFile(zip),stop=new AbortController();stop.abort();
  await assert.rejects(openProject(zip,f.session,stop.signal),/abort/i);
  assert.deepEqual(await fs.readFile(zip),before);
  assert.deepEqual(await fs.readdir(f.session),[]);
});

test('ZIP links are rejected rather than restored as filesystem references', async t => {
  const f = await fixture(t);
  const {ZipFile}=require('yazl');
  const {pipeline}=require('node:stream/promises');
  const {createWriteStream}=require('node:fs');
  const zip = new ZipFile(), file = path.join(f.temp,'link.zip');
  zip.addBuffer(Buffer.from('/outside'),'root/link',{mode:0o120777});
  const writing=pipeline(zip.outputStream,createWriteStream(file));zip.end();await writing;
  await assert.rejects(extractZip(file,f.session),/链接/);
  assert.deepEqual(await fs.readdir(f.session),[]);
});
