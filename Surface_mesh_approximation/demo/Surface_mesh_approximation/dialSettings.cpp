#include "dialSettings.h"

SettingsDialog::SettingsDialog(QWidget *parent)
  : QDialog(parent)
{
  setupUi(this);

  loadDefaults();
  loadFromSettings();
}

void SettingsDialog::loadDefaults()
{
}

void SettingsDialog::loadFromSettings()
{
  QSettings settings("settings.ini", QSettings::IniFormat);
  settings.beginGroup("PkgTSMA");

  if (settings.contains("method_random"))
    method_random->setChecked(settings.value("method_random").toBool());
  if (settings.contains("method_incremental"))
    method_incremental->setChecked(settings.value("method_incremental").toBool());
  if (settings.contains("method_hierarchical"))
    method_hierarchical->setChecked(settings.value("method_hierarchical").toBool());

  if (settings.contains("cb_nb_proxies"))
    cb_nb_proxies->setChecked(settings.value("cb_nb_proxies").toBool());
  if (settings.contains("nb_proxies"))
    nb_proxies->setValue(settings.value("nb_proxies").toInt());
  nb_proxies->setEnabled(cb_nb_proxies->isChecked());
  if (settings.contains("cb_error_drop"))
    cb_error_drop->setChecked(settings.value("cb_error_drop").toBool());
  if (settings.contains("error_drop"))
    error_drop->setValue(settings.value("error_drop").toDouble());
  error_drop->setEnabled(cb_error_drop->isChecked());

  if (settings.contains("nb_relaxations"))
    nb_relaxations->setValue(settings.value("nb_relaxations").toInt());
  if (settings.contains("nb_iterations"))
    nb_iterations->setValue(settings.value("nb_iterations").toInt());

  if (settings.contains("is_relative_to_chord"))
    is_relative_to_chord->setChecked(settings.value("is_relative_to_chord").toBool());
  if (settings.contains("with_dihedral_angle"))
    with_dihedral_angle->setChecked(settings.value("with_dihedral_angle").toBool());
  if (settings.contains("if_optimize_anchor_location"))
    if_optimize_anchor_location->setChecked(settings.value("if_optimize_anchor_location").toBool());
  if (settings.contains("pca_plane"))
    pca_plane->setChecked(settings.value("pca_plane").toBool());
  if (settings.contains("chord_error"))
    chord_error->setValue(settings.value("chord_error").toDouble());

  if (settings.contains("split_proxy_idx"))
    split_proxy_idx->setValue(settings.value("split_proxy_idx").toInt());
  if (settings.contains("split_nb_sections"))
    split_nb_sections->setValue(settings.value("split_nb_sections").toInt());
  if (settings.contains("split_nb_relaxations"))
    split_nb_relaxations->setValue(settings.value("split_nb_relaxations").toInt());

  settings.endGroup();
}

void SettingsDialog::saveToSettings()
{
  QSettings settings("settings.ini", QSettings::IniFormat);
  settings.beginGroup("PkgTSMA");

  settings.setValue("method_random", method_random->isChecked());
  settings.setValue("method_incremental", method_incremental->isChecked());
  settings.setValue("method_hierarchical", method_hierarchical->isChecked());

  settings.setValue("cb_nb_proxies", cb_nb_proxies->isChecked());
  settings.setValue("nb_proxies", nb_proxies->value());
  settings.setValue("cb_error_drop", cb_error_drop->isChecked());
  settings.setValue("error_drop", error_drop->value());

  settings.setValue("nb_relaxations", nb_relaxations->value());
  settings.setValue("nb_iterations", nb_iterations->value());

  settings.setValue("is_relative_to_chord", is_relative_to_chord->isChecked());
  settings.setValue("with_dihedral_angle", with_dihedral_angle->isChecked());
  settings.setValue("if_optimize_anchor_location", if_optimize_anchor_location->isChecked());
  settings.setValue("pca_plane", pca_plane->isChecked());
  settings.setValue("chord_error", chord_error->value());

  settings.setValue("split_proxy_idx", split_proxy_idx->value());
  settings.setValue("split_nb_sections", split_nb_sections->value());
  settings.setValue("split_nb_relaxations", split_nb_relaxations->value());

  settings.endGroup();
}

void SettingsDialog::accept()
{
  saveToSettings();

  QDialog::accept();
}
