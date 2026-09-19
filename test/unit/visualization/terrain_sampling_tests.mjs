// Run: node test/unit/visualization/terrain_sampling_tests.mjs /path/to/fixtures.json
// Fixtures are emitted by the actual-package Julia terrain export tests.
import assert from 'node:assert/strict';
import { readFile,writeFile,mkdtemp,rm } from 'node:fs/promises';
import { dirname,resolve,join } from 'node:path';
import { tmpdir } from 'node:os';
import { fileURLToPath,pathToFileURL } from 'node:url';
const repo=resolve(dirname(fileURLToPath(import.meta.url)),'../../..');
const fixture=process.argv[2]; if(!fixture)throw new Error('pass the Julia-generated fixtures.json path');
const cases=JSON.parse(await readFile(fixture,'utf8'));
assert.ok(Array.isArray(cases) && cases.length>0,'Julia terrain fixtures must be nonempty');
const combined=JSON.parse(await readFile(join(dirname(fixture),'combined-export.json'),'utf8'));
const dir=await mkdtemp(join(tmpdir(),'spaceagora-terrain-js-'));
try {
 await writeFile(join(dir,'package.json'),'{"type":"module"}');
 for(const [name,path] of [['three.js','viewer/vendor/three.module.js'],['data.js','viewer/src/data.js'],['terrain.js','viewer/src/terrain.js'],['globe.js','viewer/src/globe.js']]) {
  let text=await readFile(join(repo,path),'utf8');
  const imports={'three':'./three.js','viewer/data.js':'./data.js','viewer/globe.js':'./globe.js'};
  text=text.replace(/(['"])(three|viewer\/data\.js|viewer\/globe\.js)\1/g,(_,q,s)=>q+imports[s]+q);
  await writeFile(join(dir,name),text);
 }
  const THREE=await import(pathToFileURL(join(dir,'three.js')));
  const {FrameData}=await import(pathToFileURL(join(dir,'data.js')));
 const {createTerrain,terrainGridSample,terrainGridHeight,terrainLongitude,decodeTerrainGrids}=await import(pathToFileURL(join(dir,'terrain.js')));
  const {createGlobe,GEOMETRY_TO_BODY,globeInteriorRadiusKm}=await import(pathToFileURL(join(dir,'globe.js')));
 // Actual geometry and material creation; texture image decoding is outside Node.
 THREE.TextureLoader.prototype.load=function(url,onLoad){assert.ok(url.startsWith('data:image/jpeg;base64,'));const t=new THREE.Texture();queueMicrotask(()=>onLoad?.(t));return t;};
 const planet={equatorial_radius_m:3396200,polar_radius_m:3376000,name:'Mars',texture:'mars',spin_rad_s:[0,0,0],rotation:{t_s:[0,10],q_pi:[[0,0,0,1],[0,0,0,1]]}};
 let queries=0;
 for(const c of cases){
  const terrain=createTerrain(c.payload,planet);
  assert.equal(terrain.hasHeights,true);assert.equal(terrain.referenceRadiusKm,c.payload.reference_radius_m/1000);
  for(const q of c.queries){assert.ok(Math.abs(terrain.heightAt(q.lat,q.lon)-q.expected)<1e-4,`${c.name}: ${q.lat},${q.lon}`);queries++;}
  for(const invalid of [NaN,Infinity,-Infinity]){assert.ok(Number.isNaN(terrain.heightAt(0,invalid)));assert.ok(Number.isNaN(terrain.heightAt(invalid,0)));}
  const lat=c.payload.site.lat_deg,lon=c.payload.site.lon_deg,h=terrain.heightAt(lat,lon);
  const craftRadiusKm=terrain.referenceRadiusKm+(h+42)*1e-3;
  assert.ok(Math.abs(1000*(craftRadiusKm-terrain.referenceRadiusKm)-h-42)<1e-6);
  if(c.payload.tiles){
   assert.ok(terrain.rootMesh && terrain.hole);
   const mesh=terrain.rootMesh,positions=mesh.geometry.attributes.position;
   assert.ok(Array.from(positions.array).every(Number.isFinite));
   const middle=8*17+8,v=new THREE.Vector3().fromBufferAttribute(positions,middle).add(mesh.position);
   const root=c.payload.tiles.root;
   const centerLat=(root.lat_min+root.lat_max)/2,centerLon=(root.lon_min+root.lon_max)/2;
   const expected=c.payload.reference_radius_m+terrain.surfaceHeight(centerLat,centerLon);
   assert.ok(Math.abs(v.length()*1000-expected)<1e-5);
   assert.ok(Math.abs(terrain.siteMarker.position.length()*1000-(c.payload.reference_radius_m+h+.5))<1e-5);
   if(c.name==='imagery'){
    assert.equal(terrain.attribution,'Synthetic test image');assert.ok(terrain.modelStatus.includes('20.0 m detail'));
   }
   const globe=createGlobe(planet,null,{hole:terrain.hole});
   const shader={uniforms:{},vertexShader:'#include <common>\n#include <begin_vertex>',fragmentShader:'#include <common>\n#include <clipping_planes_fragment>'};
   globe.mesh.material.onBeforeCompile(shader);
   assert.ok(shader.fragmentShader.includes('lonDeg = mod(lonDeg, 360.0)'));
  }else{assert.equal(terrain.levels.length,0);assert.equal(terrain.hole,null);assert.equal(terrain.modelStatus,'height grids; no imagery mesh');}
 }
 console.log(`PASS Julia -> JavaScript: ${cases.length} payloads, ${queries} height queries, explicit-radius radar formula and real Three geometry`);
 const g=decodeTerrainGrids(cases[0].payload)[0];
 for(const val of [NaN,Infinity,-Infinity]){assert.ok(Number.isNaN(terrainGridSample(g,0,val)));assert.ok(Number.isNaN(terrainGridHeight(g,0,val)));assert.ok(Number.isNaN(terrainLongitude(g,val)));}
 for(const [rows,cols,data,expected] of [[1,3,[1,2,3],2],[3,1,[1,2,3],2],[1,1,[7],7]]){
  const grid=decodeTerrainGrids({grids:[{rows,cols,lat_min:0,lat_max:3,lon_min:0,lon_max:3,heights:data}]})[0];
  assert.equal(terrainGridSample(grid,1.5,1.5),expected);
  assert.equal(terrainGridSample(grid,100,-100),data[0]);
 }
 const cloned=()=>structuredClone(cases[0].payload);
 for(const value of [0,-1,1.5,true,2049]){const p=cloned();p.grids[0].rows=value;assert.throws(()=>createTerrain(p,planet));}
 for(const mutate of [p=>p.reference_radius_m=0,p=>p.reference_radius_m=Infinity,p=>p.grids[0].heights='',p=>p.grids[0].heights=[NaN],p=>p.grids[0].lon_max=p.grids[0].lon_min+360,p=>p.grids[0].lat_max=91,p=>p.grids[0].lat_min=true,p=>p.site.lon_deg=Infinity,p=>p.fallback_height_m='0']){const p=cloned();mutate(p);assert.throws(()=>createTerrain(p,planet));}
 const imageCase=cases.find(c=>c.name==='imagery');
 for(const mutate of [p=>p.tiles.nodes=[],p=>p.tiles.nodes.push({...p.tiles.nodes[0]}),p=>p.tiles.nodes[0].x=1,p=>p.tiles.nodes[0].level=21,p=>p.tiles.root.lon_max=Infinity,p=>p.tiles.max_level=1]){const p=structuredClone(imageCase.payload);mutate(p);assert.throws(()=>createTerrain(p,planet));}
 const source=await readFile(join(repo,'viewer/src/main.js'),'utf8');
 assert.ok(source.includes('if (terrain.hasHeights)'));
 assert.ok(source.includes('1000 * (r - terrain.referenceRadiusKm) - h'));
 const dust=await readFile(join(repo,'viewer/src/dust.js'),'utf8');
 assert.ok(dust.includes('referenceRadiusKm + h * DUST_M_TO_KM'));
  console.log('PASS singleton axes, malformed payload rejection, finite wrapping, imagery-free queries, main/dust reference-radius plumbing');
  // These blocks came through one actual export_visualization call, not a
  // JavaScript reconstruction of what the Julia exporter is expected to emit.
  const frames=new FrameData(combined.frames);
  assert.equal(frames.sats,2);assert.equal(frames.count,3);
  assert.deepEqual(Array.from(frames.t),[0,2,4]);
  assert.deepEqual(combined.scene.spacecraft.map(s=>s.id),[11,22]);
  assert.equal(frames.channels.length,1);assert.equal(frames.channels[0].name,'metric');
  assert.equal(frames.channelAt(0,0,0),1);assert.equal(frames.channelAt(0,4,0),5);
  assert.ok(Number.isNaN(frames.channelAt(0,2,0)));
  assert.ok(Number.isNaN(frames.channelAt(0,1,0)));assert.ok(Number.isNaN(frames.channelAt(0,3,0)));
  assert.equal(frames.channelAt(0,1,1),20);assert.equal(frames.channelAt(0,3,1),40);
  const combinedTerrain=createTerrain(combined.terrain,combined.scene.planet);
  assert.ok(combinedTerrain.hasHeights && combinedTerrain.rootMesh);
  const imageTerrain=createTerrain(imageCase.payload,combined.scene.planet);
  assert.equal(combinedTerrain.heightAt(.5,360.5),imageTerrain.heightAt(.5,360.5));
  assert.equal(combined.models['22'].url_from,'11');assert.equal(combined.models['22'].url,undefined);
  assert.ok(combined.models['11'].url.startsWith('data:'));
  assert.deepEqual([combined.models['11'].scale,combined.models['22'].scale],[1,2]);
  assert.deepEqual(combined.models['11'].rotation_deg,[0,0,0]);
  assert.deepEqual(combined.models['22'].rotation_deg,[0,0,90]);
  console.log('PASS actual combined terrain/channel/shared-model export with decimation, gaps and separate transforms');

  // A valid deep DEM remains visible from below the reference sphere. Its
  // actual bounding sphere is in the camera frustum: an assumed shallow
  // planetary floor must not incorrectly reject the entire terrain patch.
  const radiusM=1737400;
  const moon={...planet,equatorial_radius_m:radiusM,polar_radius_m:radiusM,name:'Moon',texture:'moon'};
  const deepBounds={lat_min:-.01,lat_max:.01,lon_min:-.01,lon_max:.01};
  const deepSpec={reference_radius_m:radiusM,fallback_height_m:0,
    site:{lat_deg:0,lon_deg:0,name:'synthetic deep terrain',height_m:-20000},
    grids:[{...deepBounds,rows:2,cols:2,heights:[-20000,-20000,-20000,-20000]}],
    tiles:{scheme:'quadtree',root:deepBounds,tile_px:1,max_level:0,
      nodes:[{level:0,x:0,y:0,m_per_px:1,url:'data:image/jpeg;base64,/9j/2Q=='}]}};
  const deepTerrain=createTerrain(deepSpec,moon);
  deepTerrain.group.updateMatrixWorld(true);
  const camera=new THREE.PerspectiveCamera(45,1,.001,100);
  camera.up.set(0,0,1);
  for(const cameraRadius of [radiusM/1000-4,radiusM/1000+10]){
    camera.position.set(cameraRadius,0,0);camera.lookAt(radiusM/1000-20,0,0);camera.updateMatrixWorld(true);
    const frustum=new THREE.Frustum().setFromProjectionMatrix(
      new THREE.Matrix4().multiplyMatrices(camera.projectionMatrix,camera.matrixWorldInverse));
    const bound=deepTerrain.rootMesh.geometry.boundingSphere.clone().applyMatrix4(deepTerrain.rootMesh.matrixWorld);
    assert.ok(frustum.intersectsSphere(bound),'deep terrain is independently inside the view frustum');
    deepTerrain.update(camera,900);
    assert.ok(deepTerrain.stats.drawn>0,`visible terrain culled from radius ${cameraRadius}km`);
    assert.equal(deepTerrain.rootMesh.visible,true);
  }

  // Raised terrain can be visible beyond the tangent plane of an interior
  // sphere. Check an unobstructed sightline and a genuinely hidden control.
  function terrainVisibility(height,longitude){
    const root={lat_min:-.05,lat_max:.05,lon_min:longitude-.05,lon_max:longitude+.05};
    const spec={...deepSpec,site:{...deepSpec.site,lon_deg:longitude,height_m:height},
      grids:[{...root,rows:2,cols:2,heights:Array(4).fill(height)}],tiles:{...deepSpec.tiles,root}};
    const terrain=createTerrain(spec,moon);terrain.group.updateMatrixWorld(true);
    const r=(radiusM+height)/1000,angle=THREE.MathUtils.degToRad(longitude);
    const target=new THREE.Vector3(r*Math.cos(angle),r*Math.sin(angle),0);
    const view=new THREE.PerspectiveCamera(45,1,.001,10000);
    view.up.set(0,0,1);view.position.set(radiusM/1000+100,0,0);view.lookAt(target);view.updateMatrixWorld(true);
    const frustum=new THREE.Frustum().setFromProjectionMatrix(
      new THREE.Matrix4().multiplyMatrices(view.projectionMatrix,view.matrixWorldInverse));
    const bound=terrain.rootMesh.geometry.boundingSphere.clone().applyMatrix4(terrain.rootMesh.matrixWorld);
    assert.ok(frustum.intersectsSphere(bound),'horizon fixture independently intersects the view frustum');
    terrain.update(view,900);terrain.group.updateMatrixWorld(true);
    const delta=target.clone().sub(view.position);
    const fraction=Math.max(0,Math.min(1,-view.position.dot(delta)/delta.lengthSq()));
    const closestRadius=view.position.clone().addScaledVector(delta,fraction).length();
    const ray=new THREE.Raycaster(view.position,delta.normalize(),0,10000);
    const meshes=terrain.group.children.filter(o=>o.isMesh && o.visible && o.name.startsWith('terrain-'));
    return {drawn:terrain.stats.drawn,hits:ray.intersectObjects(meshes,false).length,closestRadius};
  }
  const raised=terrainVisibility(100000,36);
  assert.ok(raised.closestRadius>radiusM/1000,'raised sightline stays outside the entire reference sphere');
  assert.ok(raised.drawn>0 && raised.hits>0,'raised visible terrain must survive horizon rejection');
  const backside=terrainVisibility(0,180);
  assert.equal(backside.drawn,0,'a patch behind the body is still rejected');

  // The shared interior bound must fit inside actual triangular globe
  // geometry, including an oblate body; analytic radii alone are insufficient.
  for(const [equatorial,polar] of [[radiusM,radiusM],[radiusM*.9,radiusM*.85]]){
    const body={...moon,equatorial_radius_m:equatorial,polar_radius_m:polar};
    const globe=createGlobe(body,null);globe.group.updateMatrixWorld(true);
    const positions=globe.mesh.geometry.attributes.position,indices=globe.mesh.geometry.index;
    const interior=globeInteriorRadiusKm(body);let triangles=0;
    assert.ok(interior>0 && Number.isFinite(interior));
    for(let i=0;i<indices.count;i+=3){
      const v=[0,1,2].map(k=>new THREE.Vector3().fromBufferAttribute(positions,indices.getX(i+k))
        .applyMatrix4(globe.mesh.matrixWorld));
      const normal=v[1].clone().sub(v[0]).cross(v[2].clone().sub(v[0]));
      if(normal.lengthSq()<1e-24)continue;
      assert.ok(interior<=Math.abs(normal.normalize().dot(v[0]))+1e-9,
        'interior sphere protrudes through a rendered globe triangle');
      triangles++;
    }
    assert.ok(triangles>0);
  }

  // Interpolation of unit vertex directions is not itself unit length. Use
  // one actual sphere triangle to distinguish a correct narrow terrain hole
  // from the old asin(interpolated.z) calculation. Node checks the injected
  // shader's dataflow and numerical contract; this is not GPU compilation.
  const sphere=createGlobe(moon,null);
  const positions=sphere.mesh.geometry.attributes.position,indices=sphere.mesh.geometry.index;
  let triangleCase;
  for(let i=0;i<indices.count;i+=3){
    const interpolated=new THREE.Vector3();
    for(let k=0;k<3;k++) interpolated.add(new THREE.Vector3()
      .fromBufferAttribute(positions,indices.getX(i+k)).normalize().applyQuaternion(GEOMETRY_TO_BODY));
    interpolated.multiplyScalar(1/3);
    const latitude=THREE.MathUtils.radToDeg(Math.asin(interpolated.z/interpolated.length()));
    const oldLatitude=THREE.MathUtils.radToDeg(Math.asin(interpolated.z));
    if(latitude>55 && latitude<65 && Math.abs(latitude-oldLatitude)>.001){
      triangleCase={interpolated,latitude,oldLatitude};break;
    }
  }
  assert.ok(triangleCase,'find a sphere triangle that discriminates fragment normalization');
  const {interpolated,latitude,oldLatitude}=triangleCase;
  const halfWidth=Math.abs(latitude-oldLatitude)/4;
  const longitude=((THREE.MathUtils.radToDeg(Math.atan2(interpolated.y,interpolated.x))%360)+360)%360;
  const hole={lat_min:latitude-halfWidth,lat_max:latitude+halfWidth,lon_min:longitude-.1,lon_max:longitude+.1};
  assert.ok(oldLatitude<hole.lat_min || oldLatitude>hole.lat_max,'old formula misses this hole');
  const withHole=createGlobe(moon,null,{hole});
  const holeShader={uniforms:{},vertexShader:'#include <common>\n#include <begin_vertex>',
    fragmentShader:'#include <common>\n#include <clipping_planes_fragment>'};
  withHole.mesh.material.onBeforeCompile(holeShader);
  const normalization=holeShader.fragmentShader.match(/\bvec3\s+(\w+)\s*=\s*normalize\(\s*vBodyDir\s*\)/);
  const latitudeInput=holeShader.fragmentShader.match(/float\s+latDeg\s*=\s*degrees\(asin\(clamp\((\w+)\.z/);
  assert.ok(latitudeInput,'recognize the fragment latitude calculation');
  const fragmentDirection=normalization && latitudeInput[1]===normalization[1]
    ? interpolated.clone().normalize() : interpolated;
  const fragmentLatitude=THREE.MathUtils.radToDeg(Math.asin(fragmentDirection.z));
  assert.ok(fragmentLatitude>=hole.lat_min && fragmentLatitude<=hole.lat_max,
    `interpolated fragment latitude ${fragmentLatitude} misses the actual triangle hole at ${latitude}`);
  assert.ok(normalization,'normalize the interpolated direction in the fragment shader');
  const longitudeInput=holeShader.fragmentShader.match(/float\s+lonDeg\s*=\s*degrees\(atan\((\w+)\.y,\s*(\w+)\.x/);
  assert.ok(longitudeInput && longitudeInput[1]===normalization[1] && longitudeInput[2]===normalization[1],
    'fragment latitude and longitude use the same normalized body direction');
  console.log('PASS deep/raised terrain visibility, hidden control, globe interior bounds and fragment-normalization contract');
} finally {await rm(dir,{recursive:true,force:true});}
