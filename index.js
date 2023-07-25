// EVENTS
window.addEventListener('DOMContentLoaded', init);

function init() {
    let onlyOne = true;

    //////////////
    // 設定用定数
    //////////////
    //// シミュレーション定数
    const areaWidth = 100; // メッシュ列数
    const areaHeight = 30; // メッシュ行数
    const dx = 1;
    const dz = 1;
    const dt = 1/1200;
    const SORIteration = 10000; // SOR法の反復回数
    const omega = 1.7; // SOR法の加速係数
    const robinA = 0.00; // ロビン境界条件の定数A
    const robinB = 0.96; // ロビン境界条件の定数B
    
    //// ナビエストークス
    const rho = 1.2928; // E-6 */ // g/mm^3 (空気の密度)
    // const Re = 1.0; // レイノルズ数 (不使用)
    const flowIO = 40000; // 外部からの流入・流出 N/mm
    
    //// パーティクル
    // 煙パーティクル
    const particleDense = 8; // 煙パーティクル密度(使ってない)
    const pAmount = 0.1; // パーティクル量
    const pTime = 1.0; // パーティクル入れる時間[s]
    const pMass = 0.4E-6; // 煙パーティクルの重さ[g]
    const pCoulomb = 10.0E-6; // 煙パーティクルの電荷量[C]
    const LvelTopVel = 2000.0; // 煙粒子が速度場から受ける力の割合
    const gX = 0; // -9.80665*1000;
    const gZ = 9.8*1000; // 重力加速度 mm/s**2
    const pSize = 5.0; // 表示サイズ
    const pColor = 0x55CCFF; // 色
    const pShift = 0.5;
    // イオンパーティクル
    const qemmitRadius = 0.25; // 荷電粒子が発生する最大半径
    const qEmmitAmount = 1200; // 荷電粒子の発生数
    // const qResist = 80;// 30.0; // 荷電粒子の空気抵抗
    const qCoulomb = 0.8E-1; // 荷電粒子の電荷量[C/m^2]
    const qMass = 4.0e-5; // 荷電粒子の質量[g/m^2]
    // const qVelToLVel = 1000; // 速度場が荷電粒子から受ける力の割合
    const qSize = 5.0; // 表示サイズ
    const qColor = 0xEE0000; // 色
    const qShift = 0.0;
    
    //// 装置の定数
    const voltage = 7000; // 印加電圧[V]
    const needleX = 0.25; // 針のX位置
    const needleHeightPix = 10; // 針と穴の距離[mm]
    const holeWidthPix = 7; // 針下の穴の半径[mm]
    const holeHeight = 0.60; // 針下の穴のZ位置
    const escHoleWidthPix = 5; // 空気穴の幅[mm] escape hole

    // 色の範囲
    const particleSize = 11; // 5; // スカラー量表示パーティクルのサイズ
    // const colors = [];
    // const color = new THREE.Color();
    // const latMinValS1 = -0.0; // 電界強度V/mm (旧圧力-5.0~5.0)
    // const latMaxValS1 = 300.0;
    // const getLatColorS1 = (val) => {return color.setHSL((1-(constrain(val,latMinValS1,latMaxValS1) - latMinValS1) / (latMaxValS1 - latMinValS1))*2/3,1,0.5);}
    // const latMinValS2 = -0.0; // 電位V
    // const latMaxValS2 = 5000.0;
    // const getLatColorS2 = (val) => {return color.setHSL((1-(constrain(val,latMinValS2,latMaxValS2) - latMinValS2) / (latMaxValS2 - latMinValS2))*2/3,1,0.5);}

    // const arrowThickness = 0.075;
    // const latMinValV1 = 10.0; // 速度 mm/s
    // const latMaxValV1 = 50.0;
    // const latV1Length = 0.01;
    // const getLatColorV1 = (val) => {return color.setHSL((1-(constrain(val,latMinValV1,latMaxValV1) - latMinValV1) / (latMaxValV1 - latMinValV1))*2/3,1,0.5);} // カラー
    // // const getLatColorV1 = (val) => {return color.setHSL(0,0,(constrain(val,latMinValV1,latMaxValV1) - latMinValV1) / (latMaxValV1 - latMinValV1));} // モノクロ
    // const latMinValV2 = 0.0; // 電界 V/mm
    // const latMaxValV2 = 300.0;
    // const latV2Length = 0.002;
    // const getLatColorV2 = (val) => {return color.setHSL((1-(constrain(val,latMinValV2,latMaxValV2) - latMinValV2) / (latMaxValV2 - latMinValV2))*2/3,1,0.5);}

    //// 色の範囲
    const colors = [];
    const latMinValS1 = -0.0; // 電界強度V/mm (旧圧力-5.0~5.0)
    const latMaxValS1 = 1000.0;
    const latMinValS2 = 0.0; // 電位V
    const latMaxValS2 = 5000.0;
    // vector
    const arrowThickness = 0.15;
    const latMinValV1 = 10.0; // 速度 mm/s
    const latMaxValV1 = 200.0;
    const latV1Length = 0.01;
    const latMinValV2 = 0.0; // 電界 V/mm
    const latMaxValV2 = 1000.0;
    const latV2Length = 0.002;
    
    // set color object's color value. referred from HSV colorbar.  arg1 .. main value | arg2 .. min | arg3 .. max
    const color = new THREE.Color();
    const getColor = (val, min, max) => {
        const range = (max - min);
        const mid = range / 2;
        const absVal = Math.abs(val-mid)/range;
        let H,S,V;
        H = S = V = 0;
        if(absVal > 0.3){
            H = (val-mid) < 0 ? 2/3 : 0;
            S = 1;
            V = 0.25 * (2-(constrain(absVal, 0.3, 0.5)-0.3)/0.2);
        } else {
            H = (1-(val - range * 0.2) / (range * 0.6))*2/3;
            S = 1;
            V = 0.5;
        }
        return color.setHSL(H,S,V);
    }
    
    // その他の定数
    let cameraDepth = 0.6; // カメラの初期深さ (大きいほど遠くからの視点) (0.5で装置の横幅に合わせる)
    const ionVel = 0;

    //// Element(ボタンや時間)の取得
    // テキストオブジェクトを取得
    var timeElement = document.querySelector("#time");
    var timeNode = document.createTextNode("");
    timeElement.appendChild(timeNode);
    timeNode.nodeValue = 0.0;
    // ボタンオブジェクトを取得
    const cbElements = {};
    const cbParams = [
        ["stop",false], // 再生＆一時停止
        ["applyVoltage",true], // 電圧印加
        ["showParticle",false], // トレーサ粒子
        ["showQ",false], // イオン粒子
        ["showWall",true], // 壁
        ["showElectrode",false], // 電極
        ["showVoltage",false], // 電圧
        ["showAbsE",false], // 電界強度
        ["showVelocity",true], // 流速
        ["showE",false], // 電界
        ["reset",false] // トレーサ粒子入れ直し
    ];
    cbParams.forEach(o=>{
        const elem = document.querySelector('#'+o[0]);
        elem.checked = true;
        cbElements[o[0]] = elem;
    });
    
    //////////////////////
    // 3D描画用オブジェクト
    //////////////////////
    //// 初期設定
    window.addEventListener('resize', onWindowResize); // EVENTS
    const renderer = new THREE.WebGLRenderer({canvas: document.querySelector('#myCanvas')}); // レンダラーを作成
    renderer.setPixelRatio(window.devicePixelRatio);
    renderer.setSize(window.innerWidth, window.innerHeight);
    const scene = new THREE.Scene(); // シーンを作成
    scene.background = new THREE.Color(0xFFFFFF); // 0xFFFFFF
    // カメラを作成
    // const camera = new THREE.PerspectiveCamera(45, areaWidth / areaHeight, 0.01, 10000); // 遠近投影カメラ
    cameraDepth *= areaWidth;
    const cameraAspect = window.innerHeight / window.innerWidth;
    const camera = new THREE.OrthographicCamera(-cameraDepth, cameraDepth, cameraAspect * cameraDepth, -cameraAspect * cameraDepth, 1, 1000); // 正投影カメラ
    const controls = new THREE.OrbitControls(camera, renderer.domElement); // カメラコントローラ適応
    // document.body.appendChild(renderer.domElement); // 謎の文(必要？)
    camera.position.set(areaWidth*dx/2, areaWidth*1.3, areaHeight*dz/2); // カメラ初期位置設定
    controls.target.set(areaWidth*dx/2, 0, areaHeight*dz/2);
    controls.update();
    const directionalLight = new THREE.DirectionalLight(0xFFFFFF); // 平行光源
    directionalLight.intensity = 1.0;
    directionalLight.position.set(1, 1, 1);
    scene.add(directionalLight);
    const ambLight = new THREE.AmbientLight(0xFFFFFF); // 環境光源
    ambLight.intensity = 1.0;
    scene.add(ambLight);

    // 円錐描画
    const geometry = new THREE.ConeGeometry(0.25, 8, 32); // (底面半径, 高さ, 底面メッシュ数)
    const material = new THREE.MeshBasicMaterial({color: 0xAAAAAA});
    const cone = new THREE.Mesh(geometry, material);
    cone.position.set(25,0.25,4); // 座標(x,y,z)
    cone.rotation.set(Math.PI/2,0,0);
    scene.add(cone);
    
    //// オブジェクトを宣言・シーンに追加
    //// Points
    // absE
    latPointGeometry1 = new THREE.BufferGeometry();
    latPointMaterial1 = new THREE.PointsMaterial({size: particleSize, vertexColors: true});
    scene.add(new THREE.Points(latPointGeometry1, latPointMaterial1));
    // Voltage
    latPointGeometry2 = new THREE.BufferGeometry();
    latPointMaterial2 = new THREE.PointsMaterial({size: particleSize, vertexColors: true});
    scene.add(new THREE.Points(latPointGeometry2, latPointMaterial2));
    //// Lines
    // // Velocity
    // latLinesGeometry1 = new THREE.BufferGeometry();
    // latLinesMaterial1 = new THREE.LineBasicMaterial({vertexColors: true});
    // scene.add(new THREE.LineSegments(latLinesGeometry1, latLinesMaterial1));
    // // E
    // latLinesGeometry2 = new THREE.BufferGeometry();
    // latLinesMaterial2 = new THREE.LineBasicMaterial({vertexColors: true});
    // scene.add(new THREE.LineSegments(latLinesGeometry2, latLinesMaterial2));
    
    //// Arrows
    // Velocity
    latArrowsGeometry1 = new THREE.BufferGeometry();
    latArrowsMaterial1 = new THREE.MeshBasicMaterial({vertexColors: true});
    scene.add(new THREE.Mesh(latArrowsGeometry1, latArrowsMaterial1));

    // E
    latArrowsGeometry2 = new THREE.BufferGeometry();
    latArrowsMaterial2 = new THREE.MeshBasicMaterial({vertexColors: true});
    scene.add(new THREE.Mesh(latArrowsGeometry2, latArrowsMaterial2));
    
    // パーティクルの設定
    const pointGeometry = new THREE.BufferGeometry();
    const pointMaterial = new THREE.PointsMaterial({color: pColor});
    pointMaterial.size = (pSize);
    scene.add(new THREE.Points(pointGeometry, pointMaterial));
    const qPointGeometry = new THREE.BufferGeometry();
    const qPointMaterial = new THREE.PointsMaterial({color: qColor});
    qPointMaterial.size = (qSize);
    scene.add(new THREE.Points(qPointGeometry, qPointMaterial));
    
    //// シミュレーション用のクラスと変数の宣言
    // 粒子クラス
    class IParticle{
        constructor(x=0, y=0, z=0, coulomb=0, mass=1, Vx=0, Vz=0){
            this.x = x;
            this.y = y;
            this.z = z;
            this.Vx = Vx;
            this.Vy = 0;
            this.Vz = Vz;
            this.mass = mass;
            this.coulomb = coulomb;
            this.existence = true;
        }
    }
    // 格子クラス
    class ILattice{
        constructor(x,y,z,barrier,metal=0){
            this.x = x;
            this.y = y;
            this.z = z;
            this.barrier = barrier;
            this.Vx = 0;
            this.Vz = 0;
            this.nextVx = 0;
            this.nextVz = 0;
            this.divergence = 0;
            this.pressure = 0;
            this.coulomb = 0;
            this.voltage = 0;
            this.Ex = 0;
            this.Ez = 0;
            this.absE = 0;
            this.cashVoltage = 0;
            this.cashEx = 0;
            this.cashEz = 0;
            this.cashAbsE = 0;
            this.metal = metal;
            this.direct = {N:0, S:0, E:0, W:0, SE:0, SW:0, NW:0, NE:0};
            this.posX = 0;
            this.posZ = 0;
            this.forceX = 0;
            this.forceZ = 0;
        }
        apV(s){
            if(s){
                this.voltage = this.cashVoltage;
                this.Ex = this.cashEx;
                this.Ez = this.cashEz;
                this.absE = this.cashAbsE;
            } else {
                this.voltage = 0;
                this.Ex = 0;
                this.Ez = 0;
                this.absE = 0;
            }
        }
    }
    ////////////////////////
    // 物理演算の定数と変数
    ////////////////////////
    // 変数
    let particle = [];
    let qParticle = [];
    const lattice = [];
    const walls = [];
    const positiveWalls = [];
    const negativeWalls = [];
    const posData = [];
    const normals = [];
    const indices = [];
    let wallMesh;
    let positiveWallMesh;
    let negativeWallMesh;
    let count;
    let errorMax;
    let u,v,t;
    const xRange = areaWidth+2;
    const zRange = areaHeight+2;
    const latticeAmount = xRange*zRange;
    let totalTime = 0;
    let centerWidth = Math.floor(areaWidth / 2);
    let centerHeight = Math.floor(areaHeight / 2);
    let needleXPix = Math.floor(areaWidth * needleX);
    let holeHeightPix = Math.floor(areaHeight * holeHeight);
    let qEmmitPos = [needleXPix,holeHeightPix-needleHeightPix];
    let invRho = 1/rho;
    const disposeArray = () => this.array = null;
    const flatParticlePos = (part) => part.map(o=>[o.x,o.y,o.z]).flat();
    const constrain = (num, low, high) => {if(num>high) return high; else if(num<low) return low; else return num}
    const posInd = (x,z) => (x<0||x>=xRange||z<0||z>=zRange) ? latticeAmount : x+z*xRange;
    const getX = (i) => i%xRange;
    const getZ = (i) => Math.floor(i/xRange);

    initial();
    tick();

    //// セットアップ関数 (電界計算もここで行う)
    function initial(){
        // 計算範囲線を作成
        const boxGeometry = new THREE.BoxGeometry(areaWidth+1.0, 0, areaHeight+1.0);
        const edges = new THREE.EdgesGeometry(boxGeometry);
        const line = new THREE.LineSegments(edges, new THREE.LineBasicMaterial({color: 0xFF00FF}));
        line.position.set(areaWidth/2+0.5, 0, areaHeight/2+0.5);
        scene.add(line);

        // 格子情報変数の生成
        for(i=0;i<latticeAmount;i++){
            const x = getX(i);
            const y = 0;
            const z = getZ(i);
            const barrier = (z==0||z==areaHeight+1)||((x==0||x==areaWidth+1)&&z>=holeHeightPix)||(z==holeHeightPix&&((x<=needleXPix-holeWidthPix||x>=needleXPix+holeWidthPix)&&x<=areaWidth-escHoleWidthPix||x==areaWidth+1)); // 壁条件を入力
            const metalP = (x==needleXPix && z==holeHeightPix-needleHeightPix); // 正電極条件を入力
            const metalN = (z==areaHeight+1)||(z==holeHeightPix&&((x<=needleXPix-holeWidthPix||x>=needleXPix+holeWidthPix)&&x<=areaWidth-escHoleWidthPix||x==areaWidth+1)); // 負電極条件を入力
            let metal = 0;
            if(metalP) metal = 1;
            else if (metalN) metal = -1;

            lattice.push(new ILattice(x, y, z, barrier, metal));
            if(barrier){
                const geometryW = new THREE.BoxGeometry(1, 0.5, 1);
                const geometryTranslated = geometryW.translate(x, 0.25, z);
                walls.push(geometryTranslated);
            }
            if(metalP){
                const geometryW = new THREE.BoxGeometry(1, 0.5, 1);
                const geometryTranslated = geometryW.translate(x, 0.25, z);
                positiveWalls.push(geometryTranslated);
            } else if(metalN){
                const geometryW = new THREE.BoxGeometry(1, 0.5, 1);
                const geometryTranslated = geometryW.translate(x, 0.25, z);
                negativeWalls.push(geometryTranslated);
            }
            lattice[i].direct = {
                N:posInd(x,z+1), S:posInd(x,z-1), E:posInd(x+1,z), W:posInd(x-1,z),
                SE:posInd(x+1,z-1), NW:posInd(x-1,z+1), NE:posInd(x+1,z+1)
            }
        }
        lattice.push(new ILattice(0,0,0,true)); // x,y,z,barrier,metal=0
        // render walls
        const wallGeometry = THREE.BufferGeometryUtils.mergeBufferGeometries(walls);
        const wallMaterial = new THREE.MeshToonMaterial({color: 0x112211});
        wallMesh = new THREE.Mesh(wallGeometry, wallMaterial);
        scene.add(wallMesh);
        const pwallGeometry = THREE.BufferGeometryUtils.mergeBufferGeometries(positiveWalls);
        const pwallMaterial = new THREE.MeshToonMaterial({color: 0x440000});
        positiveWallMesh = new THREE.Mesh(pwallGeometry, pwallMaterial);
        scene.add(positiveWallMesh);
        const nwallGeometry = THREE.BufferGeometryUtils.mergeBufferGeometries(negativeWalls);
        const nwallMaterial = new THREE.MeshToonMaterial({color: 0x000044});
        negativeWallMesh = new THREE.Mesh(nwallGeometry, nwallMaterial);
        scene.add(negativeWallMesh);
        reset();

        // 頂点インデックスの設定(矢印描画用)
        for(i=0;i<latticeAmount;i++){
            indices.push(i*3,i*3+1,i*3+2);
            normals.push(0,1,0,0,1,0,0,1,0);
        }
        latArrowsGeometry1.setIndex(indices);
        indices.length = 0;
        latArrowsGeometry1.setAttribute('normal', new THREE.Float32BufferAttribute(normals, 3).onUpload(disposeArray));
        normals.length=0;
        for(i=0;i<latticeAmount;i++){
            indices.push(i*3,i*3+1,i*3+2);
            normals.push(0,1,0,0,1,0,0,1,0);
        }
        latArrowsGeometry2.setIndex(indices);
        indices.length = 0;
        latArrowsGeometry2.setAttribute('normal', new THREE.Float32BufferAttribute(normals, 3).onUpload(disposeArray));
        normals.length=0;

        //// 電界計算 (電界計算は最初に一回だけ行う)
        // 電位の計算(反復法)
        for(t=0;t<SORIteration;t++){
            for(i=0;i<latticeAmount;i++){
                if(lattice[i].metal == 0){
                    const nextLattice = new Array(4);
                    nextLattice[0] = lattice[i].direct["E"] == latticeAmount ? robinB*lattice[i].cashVoltage : lattice[lattice[i].direct["E"]].cashVoltage;
                    nextLattice[1] = lattice[i].direct["W"] == latticeAmount ? robinB*lattice[i].cashVoltage : lattice[lattice[i].direct["W"]].cashVoltage;
                    nextLattice[2] = lattice[i].direct["N"] == latticeAmount ? robinB*lattice[i].cashVoltage : lattice[lattice[i].direct["N"]].cashVoltage;
                    nextLattice[3] = lattice[i].direct["S"] == latticeAmount ? robinB*lattice[i].cashVoltage : lattice[lattice[i].direct["S"]].cashVoltage;
                    lattice[i].cashVoltage = (1-omega) * lattice[i].cashVoltage + omega/4 * (nextLattice[0]+nextLattice[1]+nextLattice[2]+nextLattice[3]+lattice[i].coulomb);
                } else {
                    lattice[i].cashVoltage = (lattice[i].metal == 1) ? voltage : 0;
                }
            }
        }
        // 電界の計算
        for(i=0;i<latticeAmount;i++){
            lattice[i].cashEx = (lattice[i].direct["W"] == latticeAmount) ? 0 : lattice[lattice[i].direct["W"]].cashVoltage - lattice[lattice[i].direct["E"]].cashVoltage;
            lattice[i].cashEz = (lattice[i].direct["S"] == latticeAmount) ? 0 : lattice[lattice[i].direct["S"]].cashVoltage - lattice[lattice[i].direct["N"]].cashVoltage;
            lattice[i].cashAbsE = (lattice[i].cashEx**2+lattice[i].cashEz**2)**0.5;
        }
    }

    //// ループ関数
    function tick(){
        debugger
        if(!cbElements["stop"].checked){
            ////// 物理演算 //////
            //// 外力項 ////
            // 荷電粒子の生成
            if(cbElements["applyVoltage"].checked){
                for(i=0;i<100;i++) if(Math.random() < qEmmitAmount * dt) {
                    const rr = qemmitRadius;
                    const rth = Math.random();
                    qParticle.push(new IParticle(dx*(qEmmitPos[0]+rr*Math.cos(Math.PI*rth)), 0, dz*(qEmmitPos[1]+rr*Math.sin(Math.PI*rth)), qCoulomb, qMass, ionVel*(Math.random()<0.5 ? -1 : 1), ionVel));
                }
            }

            // 荷電粒子の運動と速度場への外力計算
            lattice.forEach(o => o.forceX = o.forceZ = 0);
            for(i = 0; i < qParticle.length; i++){
                const X = qParticle[i].x+qShift; // particle x+shift
                const x = Math.floor(X/dx); // lattice xIndex
                const Z = qParticle[i].z+qShift; // particle z+shift
                const z = Math.floor(Z/dz); // lattice zIndex
                const posind = posInd(x,z);
                if(posind==latticeAmount || lattice[posind].metal == -1 || Math.random() < 0.05){
                    qParticle[i].existence = false;
                } else {
                    if(cbElements["applyVoltage"].checked){
                        // const coefficient = qVelToLVel * dt;
                        const sz = Math.floor((Z-qShift) / dz); // slide z
                        let rijX = X - x;
                        let rijZ = (Z-qShift) - sz;
                        let lat = lattice[posInd(x,sz)];
                        const speedX = (((1-rijX)*lat.Ex + rijX*lattice[lat.direct["E"]].Ex)*(1-rijZ) + ((1-rijX)*lattice[lat.direct["N"]].Ex + rijX*lattice[lat.direct["NE"]].Ex)*(rijZ))*dt;
                        qParticle[i].Vx += (speedX * qParticle[i].coulomb / qParticle[i].mass) * dt; // - qResist*qParticle[i].Vx
    
                        const sx = Math.floor((X-qShift) / dx); // slide x
                        rijX = (X-qShift) - sx;
                        rijZ = Z - z;
                        lat = lattice[posInd(sx,z)];
                        const speedZ = (((1-rijX)*lat.Ez + rijX*lattice[lat.direct["E"]].Ez)*(1-rijZ) + ((1-rijX)*lattice[lat.direct["N"]].Ez + rijX*lattice[lat.direct["NE"]].Ez)*(rijZ))*dt;
                        qParticle[i].Vz += (speedZ * qParticle[i].coulomb / qParticle[i].mass) * dt; // - qResist*qParticle[i].Vz
                        if(!lattice[posind].barrier){
                            // 電界比例で外力決定
                            const a = 1.0;
                            lat.forceX += qCoulomb * (1-rijX) * (1-rijZ) * lat.Ex * a;
                            lattice[lat.direct["E"]].forceX += qCoulomb * rijX * (1-rijZ) * lattice[lat.direct["E"]].Ex * a;
                            lattice[lat.direct["N"]].forceX += qCoulomb * (1-rijX) * rijZ * lattice[lat.direct["N"]].Ex * a;
                            lattice[lat.direct["NE"]].forceX += qCoulomb * rijX * rijZ * lattice[lat.direct["NE"]].Ex * a;
                            
                            lat.forceZ += qCoulomb * (1-rijX) * (1-rijZ) * lat.Ez * a;
                            lattice[lat.direct["E"]].forceZ += qCoulomb * rijX * (1-rijZ) * lattice[lat.direct["E"]].Ez * a;
                            lattice[lat.direct["N"]].forceZ += qCoulomb * (1-rijX) * rijZ * lattice[lat.direct["N"]].Ez * a;
                            lattice[lat.direct["NE"]].forceZ += qCoulomb * rijX * rijZ * lattice[lat.direct["NE"]].Ez * a;

                            // 速度比例で外力決定
                            // const forceX = coefficient * qParticle[i].Vx;
                            // lat.forceX += forceX * (1-rijX) * (1-rijZ);
                            // lattice[lat.direct["E"]].forceX += forceX * rijX * (1-rijZ);
                            // lattice[lat.direct["N"]].forceX += forceX * (1-rijX) * rijZ;
                            // lattice[lat.direct["NE"]].forceX += forceX * rijX * rijZ;
    
                            // const forceZ = coefficient * qParticle[i].Vz;
                            // lat.forceZ += forceZ * (1-rijX) * (1-rijZ);
                            // lattice[lat.direct["E"]].forceZ += forceZ * rijX * (1-rijZ);
                            // lattice[lat.direct["N"]].forceZ += forceZ * (1-rijX) * rijZ;
                            // lattice[lat.direct["NE"]].forceZ += forceZ * rijX * rijZ;
                        }
                    }
                    qParticle[i].x += qParticle[i].Vx * dt;
                    qParticle[i].z += qParticle[i].Vz * dt;
                }
            }
            qParticle = qParticle.filter(p=>p.existence);

            //// 発散の計算 ////
            for(i=0;i<latticeAmount;i++){
                if(!lattice[i].barrier)
                    lattice[i].divergence = (-lattice[i].Vx-lattice[i].Vz+lattice[lattice[i].direct["E"]].Vx+lattice[lattice[i].direct["N"]].Vz) / dt;
            }

            //// 圧力の計算(反復法) ////
            errorMax = Infinity;
            count=0;
            while(errorMax>1&&count<1000){
                errorMax = 0;
                for(i=0;i<latticeAmount;i++){
                    const nextLattice = new Array(4);
                    if(!lattice[i].barrier){
                        nextLattice[0] = lattice[lattice[i].direct["E"]].barrier ? lattice[i].pressure : lattice[lattice[i].direct["E"]].pressure;
                        nextLattice[1] = lattice[lattice[i].direct["W"]].barrier ? lattice[i].pressure : lattice[lattice[i].direct["W"]].pressure;
                        if(getX(i)==0) nextLattice[1] = flowIO;
                        else if(getX(i)==areaWidth+1) nextLattice[0] = -flowIO;
                        nextLattice[2] = lattice[lattice[i].direct["N"]].barrier ? lattice[i].pressure : lattice[lattice[i].direct["N"]].pressure;
                        nextLattice[3] = lattice[lattice[i].direct["S"]].barrier ? lattice[i].pressure : lattice[lattice[i].direct["S"]].pressure;
                        const newPressure = (1-omega) * lattice[i].pressure + omega/4 * (nextLattice[0]+nextLattice[1]+nextLattice[2]+nextLattice[3]-rho*lattice[i].divergence);
                        errorMax = Math.max(errorMax, Math.abs(newPressure - lattice[i].pressure));
                        lattice[i].pressure = newPressure;
                    }
                }
                count++;
            }

            //// 移流項 ////
            for(i=0;i<latticeAmount;i++){
                const lat = lattice[i];
                const d = [[0,""],[0,""]];
                // x軸方向の流速の移流
                u = lat.Vx;
                v = (lattice[lat.direct["W"]].Vz + lat.Vz + lattice[lat.direct["NW"]].Vz + lattice[lat.direct["N"]].Vz)/4;
                d[0] = (u < 0) ? [-1,"E"] : [1,"W"];
                d[1] = (v < 0) ? [-1,"N"] : [1,"S"];
                lat.nextVx = -(d[0][0]*u*(lat.Vx - lattice[lat.direct[d[0][1]]].Vx) + d[1][0]*v*(lat.Vx - lattice[lat.direct[d[1][1]]].Vx))/dx;
                // z軸方向の流速の移流
                u = (lattice[lat.direct["E"]].Vx + lat.Vx + lattice[lat.direct["SE"]].Vx + lattice[lat.direct["S"]].Vx)/4;
                v = lat.Vz;
                d[0] = (u < 0) ? [-1,"E"] : [1,"W"];
                d[1] = (v < 0) ? [-1,"N"] : [1,"S"];
                lat.nextVz = -(d[0][0]*u*(lat.Vz - lattice[lat.direct[d[0][1]]].Vz) + d[1][0]*v*(lat.Vz - lattice[lat.direct[d[1][1]]].Vz))/dz;
            }

            // const aa = 50;
            // //// 電界を直接外力にする ////
            // lattice.forEach(o => {
            //     o.forceX = o.Ex * aa;
            //     o.forceZ = o.Ez * aa;
            // });

            // //// 速度の更新 ////
            // for(i=0;i<latticeAmount;i++){
            //     const lat = lattice[i];
            //     // x軸方向の流速の更新
            //     if(lattice[lat.direct["W"]].barrier||lat.barrier||lattice[lat.direct["S"]].barrier||lattice[lat.direct["N"]].barrier)
            //         lat.Vx = 0;
            //     else
            //         lat.Vx += (lat.nextVx - (lat.pressure - lattice[lat.direct["W"]].pressure) + lat.forceX) * dt;
            //     // z軸方向の流速の更新
            //     if(lattice[lat.direct["S"]].barrier||lat.barrier||lattice[lat.direct["W"]].barrier||lattice[lat.direct["E"]].barrier)
            //         lat.Vz = 0;
            //     else
            //         lat.Vz += (lat.nextVz - (lat.pressure - lattice[lat.direct["S"]].pressure) + lat.forceZ) * dt;
            // }

            //// 速度の更新(圧力項を中心差分で) ////
            for(i=0;i<latticeAmount;i++){
                const lat = lattice[i];
                // x軸方向の流速の更新
                if(lattice[lat.direct["W"]].barrier||lat.barrier||lattice[lat.direct["S"]].barrier||lattice[lat.direct["N"]].barrier)
                    lat.Vx = 0;
                else
                    lat.Vx += (lat.nextVx - (lat.pressure - lattice[lat.direct["W"]].pressure) + (lat.forceX + lattice[lat.direct["W"]].forceX)/2) * dt;
                // z軸方向の流速の更新
                if(lattice[lat.direct["S"]].barrier||lat.barrier||lattice[lat.direct["W"]].barrier||lattice[lat.direct["E"]].barrier)
                    lat.Vz = 0;
                else
                    lat.Vz += (lat.nextVz - (lat.pressure - lattice[lat.direct["S"]].pressure) + (lat.forceZ + lattice[lat.direct["S"]].forceZ)/2) * dt;
            }


            totalTime += dt;
            timeNode.nodeValue = totalTime.toFixed(2);
        }

        //// 描画
        // 微粒子の生成
        if(totalTime<pTime){
            for(i=0;i<10;i++){
                if(Math.random()<pAmount){
                    particle.push(new IParticle(dx*Math.random()+0.5, 0, dz*Math.random()*(holeHeightPix-1)+0.5, pCoulomb, pMass, 75));
                }
            }
        }

        // 粒子の運動計算と描画
        for(i = 0; i < particle.length; i++){
            if(!cbElements["stop"].checked){
                const X = particle[i].x+pShift; // particle x+shift
                const x = Math.floor(X / dx); // lattice xIndex
                const Z = particle[i].z+pShift; // particle z+shift
                const z = Math.floor(Z / dz); // lattice zIndex
                const posind = posInd(x,z);
                if(Math.abs(particle[i].x-areaWidth/2)>(areaWidth+1)/2+0.5 || Math.abs(particle[i].z-areaHeight/2)>(areaHeight+1)/2+0.5) particle[i].existence = false; // 範囲外に出たら消す
                if(lattice[posind].barrier){
                    // particle[i].existence = false; // 壁に当たったら消す
                } else {
                    const sz = Math.floor((Z-pShift) / dz); // slide z
                    let rijX = X - x;
                    let rijZ = Z-pShift - sz;
                    let lat = lattice[posInd(x,sz)];
                    const speedX = (((1-rijX)*lat.Vx + rijX*lattice[lat.direct["E"]].Vx)*(1-rijZ) + ((1-rijX)*lattice[lat.direct["N"]].Vx + rijX*lattice[lat.direct["NE"]].Vx)*(rijZ));
                    const EX = lat.Ex; //(((1-rijX)*lat.Ex + rijX*lattice[lat.direct["E"]].Ex)*(1-rijZ) + ((1-rijX)*lattice[lat.direct["N"]].Ex + rijX*lattice[lat.direct["NE"]].Ex)*(rijZ));

                    const sx = Math.floor((X-pShift) / dx); // slide x
                    rijX = X-pShift - sx;
                    rijZ = Z - z;
                    lat = lattice[posInd(sx,z)];
                    const speedZ = (((1-rijX)*lat.Vz + rijX*lattice[lat.direct["E"]].Vz)*(1-rijZ) + ((1-rijX)*lattice[lat.direct["N"]].Vz + rijX*lattice[lat.direct["NE"]].Vz)*(rijZ));
                    const EZ = lat.Ez; //(((1-rijX)*lat.Ez + rijX*lattice[lat.direct["E"]].Ez)*(1-rijZ) + ((1-rijX)*lattice[lat.direct["N"]].Ez + rijX*lattice[lat.direct["NE"]].Ez)*(rijZ));
                    
                    particle[i].Vx += (gX + LvelTopVel*(speedX-particle[i].Vx) + particle[i].coulomb * EX / particle[i].mass) * dt;
                    particle[i].Vz += (gZ + LvelTopVel*(speedZ-particle[i].Vz) + particle[i].coulomb * EZ / particle[i].mass) * dt;
                    particle[i].x += particle[i].Vx * dt;
                    particle[i].z += particle[i].Vz * dt;
                }
            }
        }
        particle = particle.filter(p => p.existence);
        pointGeometry.setAttribute('position', new THREE.Float32BufferAttribute(cbElements["showParticle"].checked ? flatParticlePos(particle) : 0, 3).onUpload(disposeArray));
        // 荷電粒子の描画
        qPointGeometry.setAttribute('position', new THREE.Float32BufferAttribute(cbElements["showQ"].checked ? flatParticlePos(qParticle) : 0, 3).onUpload(disposeArray));

        // 壁の可視化
        wallMesh.visible = cbElements["showWall"].checked;
        positiveWallMesh.visible = cbElements["showElectrode"].checked;
        negativeWallMesh.visible = cbElements["showElectrode"].checked;

        // 電界強度の可視化 ////// (旧圧力情報の可視化)
        if(cbElements["showAbsE"].checked){
            for(i=0;i<latticeAmount;i++){
                const val = lattice[i].absE;
                const x = getX(i);
                const z = getZ(i);
                posData.push(x, 0, z); // val / (latMaxValS1 - latMinValS1) * 10
                getColor(val, latMinValS1, latMaxValS1);
                colors.push(color.r,color.g,color.b);
            }
            // points position
            latPointGeometry1.setAttribute('position', new THREE.Float32BufferAttribute(posData, 3).onUpload(disposeArray));
            posData.length = 0;
            // points color
            latPointGeometry1.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3).onUpload(disposeArray));
            colors.length = 0;
        } else latPointGeometry1.setAttribute('position', new THREE.Float32BufferAttribute([],3));

        // 電位情報の可視化
        if(cbElements["showVoltage"].checked){
            for(i=0;i<latticeAmount;i++){
                const x = getX(i);
                const z = getZ(i);
                posData.push(x, 0, z); // lattice[i].voltage / (latMaxValS2 - latMinValS2) * 10
                getColor(lattice[i].voltage, latMinValS2, latMaxValS2);
                colors.push(color.r,color.g,color.b);
            }
            // points position
            latPointGeometry2.setAttribute('position', new THREE.Float32BufferAttribute(posData, 3).onUpload(disposeArray));
            posData.length = 0;
            // points color
            latPointGeometry2.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3).onUpload(disposeArray));
            colors.length = 0;
        } else latPointGeometry2.setAttribute('position', new THREE.Float32BufferAttribute([],3));
        
        // 速度情報の可視化 (xz合成ベクトル)
        // if(cbElements["showVelocity"].checked){
        //     for(i=0;i<latticeAmount;i++){
        //         const x = getX(i);
        //         const z = getZ(i);
        //         const x0 = (x - 0.5);
        //         const z0 = (z - 0.5);
        //         const dx = lattice[i].Vx*latV1Length;
        //         const dz = lattice[i].Vz*latV1Length;
        //         posData.push(x0, 0, z0, x0+dx, 0, z0+dz);
        //         getColor((lattice[i].Vx**2+lattice[i].Vz**2)**0.5, latMinValV1, latMaxValV1);
        //         colors.push(color.r,color.g,color.b);
        //         colors.push(color.r,color.g,color.b);
        //     }
        //     // lines position
        //     latLinesGeometry1.setAttribute("position",  new THREE.Float32BufferAttribute(posData, 3).onUpload(disposeArray));
        //     posData.length = 0;
        //     // lines color
        //     latLinesGeometry1.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3).onUpload(disposeArray));
        //     colors.length = 0;
        // } else latLinesGeometry1.setAttribute('position', new THREE.Float32BufferAttribute([],3));
        if(cbElements["showVelocity"].checked){
            lattice.forEach(o=>{
                const val = (o.Vx**2+o.Vz**2)**0.5;
                const val1 = o.Vx;
                const val2 = o.Vz;
                posData.push(o.x+val1 * latV1Length,0.2,o.z+val2 * latV1Length);
                if(val!=0){
                    posData.push(o.x+val2/val*arrowThickness,0.2,o.z-val1/val*arrowThickness);
                    posData.push(o.x-val2/val*arrowThickness,0.2,o.z+val1/val*arrowThickness);
                } else posData.push(o.x,0.2,o.z,o.x,0.2,o.z);
                getColor(val, latMinValV1, latMaxValV1);
                colors.push(color.r,color.g,color.b);
                colors.push(color.r,color.g,color.b);
                colors.push(0.9,0.9,0.9);
            });
            latArrowsGeometry1.setAttribute('position', new THREE.Float32BufferAttribute(posData, 3).onUpload(disposeArray));
            posData.length = 0;
            latArrowsGeometry1.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3).onUpload(disposeArray));
            colors.length = 0;
        } else latArrowsGeometry1.setAttribute('position', new THREE.Float32BufferAttribute([],3));

        // 電界情報の可視化 (xz合成ベクトル)
        // if(cbElements["showE"].checked){
        //     for(i=0;i<latticeAmount;i++){
        //         const x = getX(i);
        //         const z = getZ(i);
        //         const x0 = (x - 0.5);
        //         const z0 = (z - 0.5);
        //         posData.push(x0, 0, z0, x0+lattice[i].Ex/latV2Length, 0, z0+lattice[i].Ez/latV2Length);
        //         getColor(lattice[i].absE, latMinValV2, latMaxValV2);
        //         colors.push(color.r,color.g,color.b);
        //         colors.push(color.r,color.g,color.b);
        //     }
        //     // lines position
        //     latLinesGeometry2.setAttribute("position",  new THREE.Float32BufferAttribute(posData, 3).onUpload(disposeArray));
        //     posData.length = 0;
        //     // lines color
        //     latLinesGeometry2.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3).onUpload(disposeArray));
        //     colors.length = 0;
        // } else latLinesGeometry2.setAttribute('position', new THREE.Float32BufferAttribute([],3));
        if(cbElements["showE"].checked){
            lattice.forEach(o=>{
                const val = o.absE;
                const val1 = o.Ex;
                const val2 = o.Ez;
                posData.push(o.x+val1 * latV2Length,0.2,o.z+val2 * latV2Length);
                if(val!=0){
                    posData.push(o.x+val2/val*arrowThickness,0.2,o.z-val1/val*arrowThickness);
                    posData.push(o.x-val2/val*arrowThickness,0.2,o.z+val1/val*arrowThickness);
                } else posData.push(o.x,0.2,o.z,o.x,0.2,o.z);
                getColor(val, latMinValV2, latMaxValV2);
                colors.push(color.r,color.g,color.b);
                colors.push(color.r,color.g,color.b);
                colors.push(0.9,0.9,0.9);
            });
            latArrowsGeometry2.setAttribute('position', new THREE.Float32BufferAttribute(posData, 3).onUpload(disposeArray));
            posData.length = 0;
            latArrowsGeometry2.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3).onUpload(disposeArray));
            colors.length = 0;
        } else latArrowsGeometry2.setAttribute('position', new THREE.Float32BufferAttribute([],3));

        if(cbElements["reset"].checked) reset();
        if(cbElements["applyVoltage"].checked) lattice.forEach(o=>o.apV(true));
        else lattice.forEach(o=>o.apV(false));

        if(onlyOne){
            // 最初のフレームでオブジェクトを全体に描画しなければ、視界に原点が無い時にそのオブジェクトが描画されない。
            pointGeometry.setAttribute('position', new THREE.Float32BufferAttribute(flatParticlePos(lattice), 3).onUpload(disposeArray));
            qPointGeometry.setAttribute('position', new THREE.Float32BufferAttribute(flatParticlePos(lattice), 3).onUpload(disposeArray));
        }
        
        // required if controls.enableDamping or controls.autoRotate are set to true
        controls.update();
        // レンダリング
        renderer.render(scene, camera);
        if(onlyOne) cbParams.forEach(o=>cbElements[o[0]].checked = o[1]);
        onlyOne = false;
        requestAnimationFrame(tick);
    }

    //// 自作関数
    function reset(){
        const particleIndexW = areaWidth * particleDense;
        const totalParticles = particleIndexW * areaHeight * particleDense;
        particle.length = 0;
        totalTime = 0;
        cbElements["reset"].checked = false;
    }

    function onWindowResize(){
        const canvasWidth = window.innerWidth;
        const canvasHeight = window.innerHeight;
        renderer.setSize(canvasWidth, canvasHeight);
        const cameraAspect = window.innerHeight / window.innerWidth;
        camera.left = -cameraDepth;
        camera.right = cameraDepth;
        camera.top = cameraAspect * cameraDepth;
        camera.bottom = -cameraAspect * cameraDepth;
        camera.updateProjectionMatrix();
    }
}