// GEE 脚本：VIIRS 叶绿素浓度与实测点位进行星地匹配
// 匹配规则：时间窗口 ±1 天；空间窗口 3x3 像元（对 chlor_a 做 3x3 邻域均值）

/**** 1) 参数区（按你的资产信息修改） ****/
var inSituFc = ee.FeatureCollection('projects/ee-jiashu1225/assets/1_polar_all');

// 实测表字段名（根据你截图）
var timeField = 'DATETIME_UTC';   // 字符串时间，如 2019-12-06T09:00:00Z
var inSituChlaField = 'chla';     // 实测叶绿素浓度字段
var cruiseField = 'CRUISE';       // 航次字段（可选，仅用于保留）

// 使用 GEE 中可直接访问的海色产品：COPERNICUS/MARINE/SATELLITE_OCEAN_COLOR/V6
// 其中包含 chlor_a 以及 VIIRS_nobs（VIIRS 观测次数）
// 通过 VIIRS_nobs > 0 约束，提取“有 VIIRS 贡献”的 chlor_a 像元。
var oceanColor = ee.ImageCollection('COPERNICUS/MARINE/SATELLITE_OCEAN_COLOR/V6')
  .select(['chlor_a', 'VIIRS_nobs']);

// 时间窗口（天）
var dayWin = 1;

/**** 2) 单点匹配函数 ****/
var matchOne = function(ft) {
  var geom = ft.geometry();
  var tObs = ee.Date(ft.get(timeField));

  // ±1 天窗口，end 需要 +1 天用于包含上边界当日
  var tStart = tObs.advance(-dayWin, 'day');
  var tEnd = tObs.advance(dayWin + 1, 'day');

  // 仅保留有 VIIRS 观测的影像
  var cands = oceanColor.filterDate(tStart, tEnd)
    .filterBounds(geom)
    .map(function(img) {
      var dtHours = img.date().difference(tObs, 'hour').abs();
      var viirsChl = img.select('chlor_a').updateMask(img.select('VIIRS_nobs').gt(0));
      return viirsChl
        .copyProperties(img, img.propertyNames())
        .set('dt_hours', dtHours);
    });

  var hasImg = cands.size().gt(0);

  // 选“时间最接近”的影像
  var best = ee.Image(ee.Algorithms.If(
    hasImg,
    cands.sort('dt_hours').first(),
    null
  ));

  // 3x3 像元均值（以目标像元为中心，半径1像元）
  var sat3x3 = ee.Algorithms.If(
    hasImg,
    best.focal_mean({radius: 1, units: 'pixels'})
      .reduceRegion({
        reducer: ee.Reducer.first(),
        geometry: geom,
        scale: 4000,
        maxPixels: 1e8
      })
      .get('chlor_a'),
    null
  );

  var satOrig = ee.Algorithms.If(
    hasImg,
    best.reduceRegion({
      reducer: ee.Reducer.first(),
      geometry: geom,
      scale: 4000,
      maxPixels: 1e8
    }).get('chlor_a'),
    null
  );

  return ft.set({
    obs_time: tObs.format("YYYY-MM-dd'T'HH:mm:ss"),
    insitu_chla: ft.get(inSituChlaField),
    viirs_chla_pixel: satOrig,
    viirs_chla_3x3_mean: sat3x3,
    viirs_img_time: ee.Algorithms.If(hasImg, best.date().format("YYYY-MM-dd'T'HH:mm:ss"), null),
    time_diff_hours: ee.Algorithms.If(hasImg, best.get('dt_hours'), null),
    matched_flag: hasImg
  });
};

/**** 3) 批量匹配 ****/
var matched = inSituFc.map(matchOne);

// 预览前10条
print('Matched preview (VIIRS)', matched.limit(10));

Map.centerObject(inSituFc, 4);
Map.addLayer(inSituFc, {color: 'yellow'}, 'In-situ points');

var demoDate = ee.Date(ee.Feature(inSituFc.first()).get(timeField));
var demoImgRaw = oceanColor
  .filterDate(demoDate.advance(-1, 'day'), demoDate.advance(2, 'day'))
  .first();
var demoImg = ee.Image(demoImgRaw)
  .select('chlor_a')
  .updateMask(ee.Image(demoImgRaw).select('VIIRS_nobs').gt(0));
Map.addLayer(demoImg, {min: 0.01, max: 10, palette: ['081d58','225ea8','1d91c0','41b6c4','a1dab4','ffffcc']}, 'VIIRS chlor_a demo');

/**** 4) 导出结果到 Google Drive ****/
Export.table.toDrive({
  collection: matched,
  description: 'polar_viirs_chla_match_pm1day_3x3',
  fileFormat: 'CSV'
});

// 可选：只导出关键字段
// Export.table.toDrive({
//   collection: matched.select([
//     cruiseField, timeField, inSituChlaField,
//     'viirs_chla_pixel', 'viirs_chla_3x3_mean',
//     'viirs_img_time', 'time_diff_hours', 'matched_flag'
//   ]),
//   description: 'polar_viirs_chla_match_keyfields_pm1day_3x3',
//   fileFormat: 'CSV'
// });
