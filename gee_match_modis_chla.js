// GEE 脚本：MODIS 叶绿素浓度与实测点位进行星地匹配
// 匹配规则：时间窗口 ±1 天；空间窗口 3x3 像元（对 chlor_a 做 3x3 邻域均值）

/**** 1) 参数区（按你的资产信息修改） ****/
var inSituFc = ee.FeatureCollection('projects/ee-jiashu1225/assets/1_polar_all');

// 实测表字段名（根据你截图）
var timeField = 'DATETIME_UTC';   // 字符串时间，如 2019-12-06T09:00:00Z
var inSituChlaField = 'chla';     // 实测叶绿素浓度字段
var cruiseField = 'CRUISE';       // 航次字段（可选，仅用于保留）

// MODIS-Aqua 海色叶绿素产品
// 常用 band: chlor_a
var modis = ee.ImageCollection('NASA/OCEANDATA/MODIS-Aqua/L3SMI')
  .select('chlor_a');

// 时间窗口（天）
var dayWin = 1;

/**** 2) 单点匹配函数 ****/
var matchOne = function(ft) {
  var geom = ft.geometry();
  var tObs = ee.Date(ft.get(timeField));

  // ±1 天窗口，end 需要 +1 天用于包含上边界当日
  var tStart = tObs.advance(-dayWin, 'day');
  var tEnd = tObs.advance(dayWin + 1, 'day');

  var cands = modis.filterDate(tStart, tEnd)
    .filterBounds(geom)
    .map(function(img) {
      var dtHours = img.date().difference(tObs, 'hour').abs();
      return img.set('dt_hours', dtHours);
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
        scale: 4638,      // MODIS L3SMI 典型分辨率 ~4.6km
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
      scale: 4638,
      maxPixels: 1e8
    }).get('chlor_a'),
    null
  );

  return ft.set({
    obs_time: tObs.format("YYYY-MM-dd'T'HH:mm:ss"),
    insitu_chla: ft.get(inSituChlaField),
    modis_chla_pixel: satOrig,         // 原像元值
    modis_chla_3x3_mean: sat3x3,       // 3x3 空间窗口均值
    modis_img_time: ee.Algorithms.If(hasImg, best.date().format("YYYY-MM-dd'T'HH:mm:ss"), null),
    time_diff_hours: ee.Algorithms.If(hasImg, best.get('dt_hours'), null),
    matched_flag: hasImg
  });
};

/**** 3) 批量匹配 ****/
var matched = inSituFc.map(matchOne);

// 预览前10条
print('Matched preview', matched.limit(10));

// 在地图上看一下点位与示例影像
Map.centerObject(inSituFc, 4);
Map.addLayer(inSituFc, {color: 'yellow'}, 'In-situ points');

var demoDate = ee.Date(ee.Feature(inSituFc.first()).get(timeField));
var demoImg = modis.filterDate(demoDate.advance(-1, 'day'), demoDate.advance(2, 'day')).first();
Map.addLayer(demoImg, {min: 0.01, max: 10, palette: ['081d58','225ea8','1d91c0','41b6c4','a1dab4','ffffcc']}, 'MODIS chlor_a demo');

/**** 4) 导出结果到 Google Drive ****/
Export.table.toDrive({
  collection: matched,
  description: 'polar_modis_chla_match_pm1day_3x3',
  fileFormat: 'CSV'
});

// 可选：只导出关键字段
// Export.table.toDrive({
//   collection: matched.select([
//     cruiseField, timeField, inSituChlaField,
//     'modis_chla_pixel', 'modis_chla_3x3_mean',
//     'modis_img_time', 'time_diff_hours', 'matched_flag'
//   ]),
//   description: 'polar_modis_chla_match_keyfields_pm1day_3x3',
//   fileFormat: 'CSV'
// });
