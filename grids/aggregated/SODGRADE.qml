<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="2.14.1-Essen" minimumScale="0" maximumScale="1e+08" hasScaleBasedVisibilityFlag="0">
  <pipe>
    <rasterrenderer opacity="1" alphaBand="-1" classificationMax="4" classificationMinMaxOrigin="User" band="1" classificationMin="0" type="singlebandpseudocolor">
      <rasterTransparency/>
      <rastershader>
        <colorrampshader colorRampType="INTERPOLATED" clip="0">
          <item alpha="255" value="0" label="No sodic properties" color="#f7f4f9"/>
          <item alpha="255" value="1" label="Some sodic properties" color="#d4b9da"/>
          <item alpha="255" value="2" label="Moderate sodic properties" color="#df65b0"/>
          <item alpha="255" value="3" label="Solonetz or Calcisols, sodic properties and/or ph >8.1" color="#ce1256"/>
          <item alpha="255" value="4" label="Solonetz, Solonchak, gypsic properties and/or ph >8.5" color="#67001f"/>
        </colorrampshader>
      </rastershader>
    </rasterrenderer>
    <brightnesscontrast brightness="0" contrast="0"/>
    <huesaturation colorizeGreen="128" colorizeOn="0" colorizeRed="255" colorizeBlue="128" grayscaleMode="0" saturation="0" colorizeStrength="100"/>
    <rasterresampler maxOversampling="2"/>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
